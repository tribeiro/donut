#! /usr/bin/env python

import os,sys
import numpy as np
import pylab as py
from astropy import units as u
from astropy.coordinates import Angle
import logging

logging.basicConfig(format='%(asctime)s[%(levelname)s]-%(name)s-(%(filename)s:%(lineno)d):: %(message)s',
                    level=logging.INFO)

log = logging.getLogger(__name__)


def main(argv):

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-f','--filename',
                      help = 'Input file with output from donut (.npy).'
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
    cat = np.load(opt.filename).T

    pix2mm = 0.01 # pixel size in um
    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commax = 8
    id_commay = 9

    ####################################################################################################################

    def fit():

        mask = np.zeros_like(x) == 0
        nreject = 0

        for iter in range(niter):
            x1,y1 = x[mask],y[mask]

            z = np.polyfit(x1,y1,1)

            fit = np.poly1d(z)

            yres = y - fit(x) # residue
            avg = np.mean(yres)
            std = np.std(yres)
            mask = np.sqrt(yres**2.) < std*1.5
            new_nreject = len(mask)-len(mask[mask])

            log.debug('Iter[%i/%i]: Avg = %f Std = %f Reject. %i'%(iter,
                                                                  niter,
                                                                  avg,
                                                                  std,
                                                                  nreject))

            if new_nreject > opt.maxreject:
                log.debug('Maximum reject reached. Breaking.')
                break
            elif new_nreject == nreject:
                log.debug('Rejected 0 data points. Breaking')
                break
            nreject = new_nreject

        return z,mask

    ####################################################################################################################

    def plot():

        py.plot(x,
                y,'o',markerfacecolor='w')
        py.plot(x[mask],
                y[mask],'bo')

        py.plot(xx,
                newFit(xx),'r-')

        ylim = py.ylim()
        xlim = py.xlim()

        if xlim[0] < root < xlim[1]:
            py.plot([root,root],ylim,'r--')
            py.ylim(ylim)

        mean = np.mean(y[mask])
        py.plot(xlim,[mean,mean],
                'b--')

        py.grid()

        return

    ####################################################################################################################

    center = [9216/2,9232/2]
    xx = (cat[0]-center[0])
    yy = (cat[1]-center[1])
    pre_mask = np.bitwise_and(np.abs(xx) > 100, np.abs(yy)> 100)

    xx = (cat[0][pre_mask]-center[0])
    yy = (cat[1][pre_mask]-center[1])

    x = np.sqrt(xx**2.+yy*2.)*pix2mm #* yy/np.abs(yy)
    y = cat[id_astigx][pre_mask]

    niter = opt.niter if opt.niter > 0 else 1

    z,mask = fit()
    root = -z[1]/z[0]
    V_angle = Angle(z[0]*u.rad)

    newFit = np.poly1d(z)
    mean = np.mean(y[mask])
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('Primary mirror astigmatism in V = %f'%(mean/10.))
    log.info('Angle V = %s degrees'%(V_angle.to_string(unit=u.degree, sep=(':', ':', ':'))))
    log.info('Center X = %.4f mm / %.2f pixels'%(root*1e-3, root/pix2mm))


    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.subplot(231)

    plot()

    ####################################################################################################################

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_astigy][pre_mask]

    z,mask = fit()
    root = -z[1]/z[0]
    U_angle = Angle(z[0]*u.rad)
    newFit = np.poly1d(z)
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    mean = np.mean(y[mask])
    log.info('Primary mirror astigmatism in  U = %f'%(mean/10.))
    log.info('Angle U = %s degrees'%(U_angle.to_string(unit=u.degree, sep=(':', ':', ':'))))
    log.info('Center Y = %.4f mm / %.2f pixels'%(root*1e-3, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.subplot(232)

    plot()

    ####################################################################################################################

    py.subplot(233)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_seeing][pre_mask]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    Seeing = z[1]
    log.info('Seeing = %.4f arcsec'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    plot()

    ####################################################################################################################

    py.subplot(234)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_commax][pre_mask]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    ShiftX = z[1]
    log.info('ShiftX = %.4f mm'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    plot()

    ####################################################################################################################

    py.subplot(235)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_commay][pre_mask]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    ShiftY = z[1]
    log.info('ShiftY = %.4f mm'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    plot()

    ####################################################################################################################

    py.subplot(236)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_focus][pre_mask]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d([z[0],0])
    log.debug('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    yfoc = y-newFit(x)
    ShiftZ = np.mean(yfoc[mask])
    log.info('ShiftZ = %.4f mm'%(ShiftZ))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            yfoc,'o',markerfacecolor='w')
    py.plot(x[mask],
            yfoc[mask],'bo')

    py.plot(xx,
            np.zeros_like(xx)+np.mean(yfoc[mask]),'r-')

    # ylim = py.ylim()
    #
    # py.plot([root,root],ylim,'r--')
    # py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    # Z4 + h [2 - cos(alpha_x) - cos(alpha_y)]
    CFP = 291.36 # Comma Free Point in mm
    zhexapod = ShiftZ #- CFP * (2. - np.cos(V_angle.rad) - np.cos(U_angle.rad))
    xhexapod = ShiftX #+ CFP * np.sin(V_angle.rad)
    yhexapod = ShiftY #+ CFP * np.sin(U_angle.rad)

    print('''Hexapod offset:
    X = %+8.4f um
    Y = %+8.4f um
    Z = %+8.4f um
    U = %s degrees
    V = %s degrees
    '''%(xhexapod/10,
         yhexapod/10,
         zhexapod/10,
         U_angle.to_string(unit=u.degree, sep=(':', ':', ':')),
         V_angle.to_string(unit=u.degree, sep=(':', ':', ':')))
             )
    py.show()

if __name__ == '__main__':

    main(sys.argv)