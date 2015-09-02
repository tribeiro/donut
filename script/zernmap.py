#! /usr/bin/env python

import os,sys
import numpy as np
from scipy.interpolate import Rbf,bisplev,bisplrep
import pylab as py
from astropy import units as u
from astropy.coordinates import Angle
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

    pix2mm = 0.01 # pixel size in um
    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commay = 8
    id_commax = 9

    def map():
        mean = np.mean(z)
        std = np.std(z)
        mask = np.bitwise_and(z < (mean+1*std),
                              z > (mean-1*std))

        new_x = np.linspace(x.min(),x.max(),101)
        new_y = np.linspace(y.min(),y.max(),101)

        Xi,Yi = np.meshgrid(new_x,new_y)

        # try:
        # rbf = Rbf(x[mask], y[mask], z[mask])#,epsilon=4)
        # Zi = rbf(Xi,Yi)
        tck,fp,ier,msg = bisplrep(x[mask], y[mask], z[mask],kx=3,ky=3,full_output=1)
        #f = interp2d(x[mask], y[mask], z[mask],kind='cubic')
        Zi = bisplev(new_x,new_y,tck)
        # except np.linalg.linalg.LinAlgError,e:
        #     Zi = np.zeros_like(Xi)
        #     log.exception(e)
        #     pass
        # except Exception,e:
        #     raise

        log.debug('Plotting...')

        # print np.min(z[mask]),np.max(z[mask])
        py.scatter(x,y,100,z,marker='p',vmin=np.min(z[mask]),vmax=np.max(z[mask]))
        # map = py.imshow(Zi,origin='lower',vmin=np.min(z[mask]),vmax=np.max(z[mask]),extent=[new_x[0],new_x[-1],new_y[0],new_y[-1]])
        py.xlim(0,9216)
        py.ylim(0,9232)
        py.colorbar()
        log.debug('Done')

    log.debug('Mapping astigmatism in 0 degrees.')

    x = cat[0]
    y = cat[1]
    z = cat[id_astigx]

    fig = py.figure(1)

    py.subplot(231)
    # py.subplot(111)
    map()
    # py.show()
    # return

    log.debug('Mapping astigmatism in 45 degrees.')
    x = cat[0]
    y = cat[1]
    z = cat[id_astigy]
    py.subplot(232)
    map()

    log.debug('Mapping total astigmatism.')
    x = cat[0]
    y = cat[1]
    z = np.sqrt(cat[id_astigy]**2 + cat[id_astigx]**2)
    index = np.argmin(z)
    log.debug('Minimum astigmatism is %f @ %f x %f'%(z[index],x[index],y[index]))
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
    z = cat[id_commax]
    py.subplot(234)
    map()

    log.debug('Mapping comma in Y.')
    z = cat[id_commay]
    py.subplot(235)
    map()

    log.debug('Mapping total comma.')
    z = np.sqrt(cat[id_commay]**2.+cat[id_commax]**2.)
    py.subplot(236)
    map()

    # log.debug('Mapping focus.')
    # z = cat[id_focus]
    # py.subplot(236)
    # map()


    py.show()

if __name__ == '__main__':
    main(sys.argv)