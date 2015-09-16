#! /usr/bin/env python

import os,sys
from collections import OrderedDict
import numpy as np
from scipy.interpolate import Rbf,bisplev,bisplrep
from scipy.optimize import leastsq
import pylab as py
from astropy import units as u
from astropy.coordinates import Angle
from astropy.io import fits
import logging

logging.basicConfig(format='%(asctime)s[%(levelname)s]-%(name)s-(%(filename)s:%(lineno)d):: %(message)s',
                    level=logging.DEBUG)

log = logging.getLogger(__name__)


def main(argv):

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-f','--filename',
                      help = 'Input file with output from donut (.npy).'
                      ,type='string')
    parser.add_option('-i','--image',
                      help = 'Fits file where the zernike coefficients where measured (to read current offsets).'
                      ,type='string')
    parser.add_option('--niter',
                      help = 'Number of iterations on the linear fitting procedure (default = 1).'
                      ,type='int',default=1)
    parser.add_option('--maxreject',
                      help = 'Maximum number of rejected points (default=3).'
                      ,type='int',default=3)
    parser.add_option('-o','--output',
                      help = 'Output file name.'
                      ,type='string')

    opt, args = parser.parse_args(argv)

    log.info('Reading input catalog: %s'%opt.filename)
    rcat = np.load(opt.filename).T

    cat = np.array([])

    for col in rcat:
        cat = np.append(cat,np.array([col.reshape(-1),]))
    cat = cat.reshape(rcat.shape[0],-1)
    fitmask = cat[0] == 1
    cat = cat[1:]

    niter = opt.niter if opt.niter > 0 else 1
    pix2mm = 0.01 # pixel size in um
    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commay = 8
    id_commax = 9

    def fitPlaneOptimize(XYZ):
        def residiuals(parameter,f,x,y):
            return [(f[i] - model(parameter,x[i],y[i])) for i in range(len(f))]


        def model(parameter, x, y):
            a, b, c = parameter
            return a*x + b*y + c

        X = XYZ[:,0]
        Y = XYZ[:,1]
        Z = XYZ[:,2]

        p0 = [1., 1.,1.] # initial guess
        result = leastsq(residiuals, p0, args=(Z,X,Y))[0]

        return result

    def astigmatism(Axis=0):

        median = np.median(z)
        std = np.std(z)
        mask = np.abs(z-median) < std*2
        # print median,std
        fitx = (x[mask]-center[0])*pix2mm
        fity = (y[mask]-center[1])*pix2mm
        plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)
        ix = Angle(plane[0]*u.rad)/10.
        iy = Angle(plane[1]*u.rad)/10.
        CFP = 291.36 # Comma Free Point in mm
        zcfp = - CFP * (2. - np.cos(ix.rad) - np.cos(iy.rad))
        xcfp = CFP * np.sin(iy.rad)
        ycfp = CFP * np.sin(ix.rad)

        res = OrderedDict([('X',0),
                           ('Y',0),
                           ('Z',0),
                           ('U',0),
                           ('V',0)])

        if Axis == 0:
            res['U'] = -ix
            res['V'] = iy
            res['X'] = xcfp
            res['Y'] = ycfp
            res['Z'] = zcfp
        else:
            res['U'] = iy
            res['V'] = ix
            res['X'] = ycfp
            res['Y'] = -xcfp
            res['Z'] = zcfp
        print "###########################"
        print "Inclination U: %s"%res['U'].to_string(unit=u.degree, sep=(':', ':', ' '))
        print "Inclination V: %s"%res['V'].to_string(unit=u.degree, sep=(':', ':', ' '))
        print 'Corrections:'
        print 'X: ',res['X']
        print 'Y: ',res['Y']
        print 'Z: ',res['Z']
        print "###########################"

        newx = np.linspace(-center[0]*pix2mm,center[0]*pix2mm,101)
        newy = np.linspace(-center[1]*pix2mm,center[1]*pix2mm,101)
        XX,YY = np.meshgrid(newx,newy)
        ZZ = plane[0]*XX + plane[1]*YY + plane[2]
        py.pcolor(XX/pix2mm+center[0],YY/pix2mm+center[1],ZZ,vmin=zmin,vmax=zmax)
        # py.colorbar()
        py.scatter(x[mask],y[mask],50,z[mask],marker='o',vmin=zmin,vmax=zmax)
        py.xlim(0,9216)
        py.ylim(0,9232)
        py.colorbar()

        return res

    ####################################################################################################################

    def comma(Axis = 0):

        mask = np.zeros_like(x) == 0
        nreject = 0
        pos =  (x[mask]-center[0])*pix2mm  if Axis == 0 else (y[mask]-center[1])*pix2mm

        mask = np.bitwise_and(pos > -10. ,
                              pos <  10.  )

        zero_reject = len(mask)-len(mask[mask])
        nreject = len(mask)-len(mask[mask])
        for iter in range(niter):
            x1,y1 = pos[mask],z[mask]

            pol = np.polyfit(x1,y1,1)

            fit = np.poly1d(pol)

            yres = z - fit(pos) # residue
            avg = np.mean(yres[mask])
            std = np.std(yres[mask])
            new_mask = np.sqrt((yres-avg)**2.) < std
            mask = np.bitwise_and(mask,new_mask)
            new_nreject = len(mask)-len(mask[mask])

            log.debug('Iter[%i/%i]: Avg = %f Std = %f Reject. %i'%(iter,
                                                                  niter,
                                                                  avg,
                                                                  std,
                                                                  new_nreject))

            if new_nreject-zero_reject > opt.maxreject:
                log.debug('Maximum reject (%i) reached (%i). Breaking.'%(opt.maxreject,new_nreject))
                break
            elif new_nreject == nreject:
                log.debug('Rejected 0 data points. Breaking')
                break
            nreject = new_nreject

        return pol,mask

    ####################################################################################################################

    def plot(mask,Axis=0):

        pos =  (x-center[0])*pix2mm  if Axis == 0 else (y-center[1])*pix2mm
        val = z

        py.plot(pos[mask],
                val[mask],'bo')
        ylim = py.ylim()
        xlim = py.xlim()
        print xlim
        py.plot(pos[np.bitwise_not(mask)],
                val[np.bitwise_not(mask)],'o',markerfacecolor='w')

        xx =  (x[mask]-center[0])*pix2mm  if Axis == 0 else (y[mask]-center[1])*pix2mm
        py.plot(xx,
                newFit(xx),'r-')
        py.plot(xx,
                newFit(xx)+np.std(val[mask]),'r:')
        py.plot(xx,
                newFit(xx)-np.std(val[mask]),'r:')

        # ylim = py.ylim()


        if xlim[0] < root < xlim[1]:
            py.plot([root,root],ylim,'r--')
        py.plot(xlim,[newFit(0.),newFit(0.)],'k-')
        py.plot([0,0],ylim,'k--')
        py.ylim(ylim)
        py.xlim(xlim)

        mean = np.mean(y[mask])
        py.plot(xlim,[mean,mean],
                'b--')

        py.grid()

        return

    ####################################################################################################################

    def map():
        mean = np.mean(z)
        std = np.std(z)
        mask = np.bitwise_and(z < (mean+3*std),
                              z > (mean-3*std))

        if len(mask[mask]) < 100:
            mask = np.zeros_like(z) == 0

        new_x = np.linspace(x.min(),x.max(),101)
        new_y = np.linspace(y.min(),y.max(),101)

        Xi,Yi = np.meshgrid(new_x,new_y)

        # try:
        # rbf = Rbf(x[mask], y[mask], z[mask])#,epsilon=4)
        # Zi = rbf(Xi,Yi)
        # tck,fp,ier,msg = bisplrep(x[mask], y[mask], z[mask],kx=3,ky=3,full_output=1)
        #f = interp2d(x[mask], y[mask], z[mask],kind='cubic')
        # Zi = bisplev(new_x,new_y,tck)
        fitx = (x[mask]-center[0])*pix2mm
        fity = (y[mask]-center[1])*pix2mm

        try:
            plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)
            newx = np.linspace(-center[0]*pix2mm,center[0]*pix2mm,101)
            newy = np.linspace(-center[1]*pix2mm,center[1]*pix2mm,101)
            XX,YY = np.meshgrid(newx,newy)
            ZZ = plane[0]*XX + plane[1]*YY + plane[2]
            py.pcolor(XX,YY,ZZ,vmin=np.min(z[mask]),vmax=np.max(z[mask]))

        except:
            pass
        # except np.linalg.linalg.LinAlgError,e:
        #     Zi = np.zeros_like(Xi)
        #     log.exception(e)
        #     pass
        # except Exception,e:
        #     raise

        log.debug('Plotting...')

        # print np.min(z[mask]),np.max(z[mask])
        # print [new_x[0],new_x[-1],new_y[0],new_y[-1]]
        # map = py.imshow(Zi,origin='lower',vmin=np.min(z[mask]),vmax=np.max(z[mask]),extent=[new_x[0],new_x[-1],new_y[0],new_y[-1]])

        py.scatter(fitx,fity,100,z[mask],marker='p',vmin=np.min(z[mask]),vmax=np.max(z[mask]))
        py.xlim(-center[0]*pix2mm,center[0]*pix2mm)
        py.ylim(-center[1]*pix2mm,center[1]*pix2mm)
        py.colorbar()
        log.debug('Done')

        return plane

    log.debug('Mapping astigmatism...')

    fig = py.figure(1)

    py.subplot(231)

    zmin,zmax = -0.05,0.05


    center = [9216/2,9232/2]
    x = cat[0]
    y = cat[1]
    z = cat[id_astigx]

    py.subplot(231)
    planeU = astigmatism(0)

    z = cat[id_astigy]

    py.subplot(232)
    planeV = astigmatism(1)

    log.debug('Mapping Seeing.')
    x = cat[0]
    y = cat[1]
    z = cat[id_seeing]

    py.subplot(233)
    map()
    # py.show()
    # return
    #
    # log.debug('Mapping Seeing.')
    # z = cat[id_seeing]
    # py.subplot(233)
    # map()

    log.debug('Mapping comma in X.')
    x = cat[0]
    y = cat[1]
    z = cat[id_commax]

    py.subplot(234)
    cmaX,mask = comma(0)
    root = -cmaX[1]/cmaX[0]
    newFitX = np.poly1d(cmaX)
    newFit = newFitX
    plot(mask,0)

    log.debug('Mapping comma in Y.')
    z = cat[id_commay]
    py.subplot(235)
    cmaY,mask = comma(1)
    root = -cmaY[1]/cmaY[0]
    newFitY = np.poly1d(cmaY)
    newFit = newFitY
    plot(mask,1)

    log.debug('Mapping focus.')
    z = cat[id_focus]
    py.subplot(236)
    focus = map()

    U = (planeU['U']+planeV['U'])/2.
    V = (planeU['V']+planeV['V'])/2.

    print '#'*39
    print '# Offset X: %+6.4f (%+6.4f/%+6.4f) #'%(-newFitX(0.)+(planeU['X']+planeV['X'])/2.,-newFitX(0.),(planeU['X']+planeV['X'])/2.)
    print '# Offset Y: %+6.4f (%+6.4f/%+6.4f) #'%(-newFitY(0.)+(planeU['Y']+planeV['Y'])/2.,-newFitY(0.),(planeU['Y']+planeV['Y'])/2.)
    print '# Offset Z: %+6.4f %s #'%(focus[2]/10.,' '*17)
    print '# Offset U: %25s #'%(U.to_string(unit=u.degree, sep=(':', ':', ' ')))
    print '# Offset V: %25s #'%(V.to_string(unit=u.degree, sep=(':', ':', ' ')))


    if opt.image is not None:
        print '#'*56
        hdr = fits.getheader(opt.image)
        print '# Offset X: %+6.4f%+6.4f = %+6.4f (%+6.4f/%+6.4f) #'%(hdr['DXHEX'],-newFitX(0.)+(planeU['X']+planeV['X'])/2.,
                                                              hdr['DXHEX']-newFitX(0.)+(planeU['X']+planeV['X'])/2.,
                                                              hdr['DXHEX']-newFitX(0.),hdr['DXHEX']+(planeU['X']+planeV['X'])/2.)
        print '# Offset Y: %+6.4f%+6.4f = %+6.4f (%+6.4f/%+6.4f) #'%(hdr['DYHEX'],-newFitY(0.)+(planeU['Y']+planeV['Y'])/2.,
                                                                     hdr['DYHEX']-newFitY(0.)+(planeU['Y']+planeV['Y'])/2.,
                                                                     hdr['DYHEX']-newFitY(0.),hdr['DYHEX']+(planeU['Y']+planeV['Y'])/2.)
        print '# Offset Z: %+6.4f%+6.4f = %+6.4f %s #'%(hdr['DZHEX'],focus[2]/10.,hdr['DZHEX']+(focus[2]/10.),' '*17)

        du = Angle(hdr['DUHEX']*u.degree)
        dv = Angle(hdr['DVHEX']*u.degree)
        corrU = du+U
        corrV = dv+V
        print '# Offset U: %s%s = %s    #'%(du.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                          U.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                          corrU.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True))
        print '# Offset V: %s%s = %s    #'%(dv.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                          V.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                          corrV.to_string(unit=u.degree, sep=':',precision=2,alwayssign=True,pad=True))
        print '#'*56
    else:
        print '#'*39
    # log.debug('Mapping focus.')
    # z = cat[id_focus]
    # py.subplot(236)
    # map()


    py.show()

if __name__ == '__main__':
    main(sys.argv)