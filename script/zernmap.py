#! /usr/bin/env python

import os,sys
import numpy as np
from scipy.interpolate import Rbf
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
    id_commax = 8
    id_commay = 9

    def map():
        new_x = np.linspace(x.min(),x.max(),2*len(x))
        new_y = np.linspace(y.min(),y.max(),2*len(y))

        Xi,Yi = np.meshgrid(new_x,new_y)

        rbf = Rbf(x, y, z, epsilon=2)
        Zi = rbf(Xi,Yi)

        py.pcolor(Xi, Yi, Zi)
        py.scatter(x,y,100,z)
        py.xlim(0,9216)
        py.ylim(0,9232)

    log.debug('Mapping astigmatism in X.')

    x = cat[0]
    y = cat[1]
    z = cat[id_astigx]

    py.subplot(231)
    map()

    log.debug('Mapping astigmatism in Y.')
    z = cat[id_astigy]
    py.subplot(232)
    map()

    log.debug('Mapping Seeing.')
    z = cat[id_seeing]
    py.subplot(233)
    map()

    log.debug('Mapping comma in X.')
    z = cat[id_commax]
    py.subplot(234)
    map()

    log.debug('Mapping comma in Y.')
    z = cat[id_commay]
    py.subplot(235)
    map()

    log.debug('Mapping focus.')
    z = cat[id_focus]
    py.subplot(236)
    map()


    py.show()

if __name__ == '__main__':
    main(sys.argv)