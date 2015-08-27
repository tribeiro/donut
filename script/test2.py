#! /usr/bin/env python

import os,sys
import numpy as np
import pylab as py
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

    pix2mm = 0.009 # pixel size in mm
    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commax = 8
    id_commay = 9

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
            mask = np.sqrt(yres**2.) < std*3
            new_nreject = len(mask)-len(mask[mask])

            log.info('Iter[%i/%i]: Avg = %f Std = %f Reject. %i'%(iter,
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

    x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_astigx]

    niter = opt.niter if opt.niter > 0 else 1

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('Angle = %.4f degrees'%(z[0]*180/np.pi))
    log.info('Center X = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.subplot(231)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    ylim = py.ylim()

    py.plot([root,root],ylim,'r--')
    py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.subplot(232)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_astigy]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('Angle = %.4f degrees'%(z[0]*180/np.pi))
    log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    ylim = py.ylim()

    py.plot([root,root],ylim,'r--')
    py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.subplot(233)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_seeing]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('Seeing = %.4f arcsec'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    # ylim = py.ylim()

    # py.plot([root,root],ylim,'r--')
    # py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.subplot(234)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_commax]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('ShiftX = %.4f mm'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    ylim = py.ylim()

    py.plot([root,root],ylim,'r--')
    py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.subplot(235)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_commay]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('ShiftY = %.4f mm'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    ylim = py.ylim()

    py.plot([root,root],ylim,'r--')
    py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.subplot(236)

    # x = np.sqrt(cat[0]**2. + cat[1]**2.)*pix2mm
    y = cat[id_focus]

    z,mask = fit()
    root = -z[1]/z[0]
    newFit = np.poly1d(z)
    log.info('Rejected %i of %i'%(len(mask)-len(mask[mask]),
                                  len(mask)))
    log.info('ShiftZ = %.4f mm'%(z[1]))
    #log.info('Center Y = %.4f mm / %.2f pixels'%(root, root/pix2mm))

    xx = np.linspace(x.min()-200*pix2mm,x.max()+200*pix2mm)

    py.plot(x,
            y,'o',markerfacecolor='w')
    py.plot(x[mask],
            y[mask],'bo')

    py.plot(xx,
            newFit(xx),'r-')

    # ylim = py.ylim()
    #
    # py.plot([root,root],ylim,'r--')
    # py.ylim(ylim)

    py.grid()

    ####################################################################################################################

    py.show()

if __name__ == '__main__':

    main(sys.argv)