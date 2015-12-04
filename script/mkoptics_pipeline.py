#!/usr/env python

'''
Pipeline for determining hexapod corrections.
'''

import sys,os
import numpy as np
from chimera.util.sextractor import SExtractor
import logging
import subprocess
from astropy import units
from astropy.coordinates import Angle
from astropy.io import fits
from donut.zernmap import ZernMap
# import zernmap

log = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s[%(levelname)s:%(threadName)s]-%(name)s-(%(filename)s:%(lineno)d):: %(message)s',
                    level=logging.DEBUG)
def main(argv):

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-i','--image',
                      help = 'Fit image file with donuts to measure.'
                      ,type='string')
    parser.add_option('-p','--parameters',
                  help = 'Base donut json parameters. Position will be overwritten.'
                  ,type='string')
    parser.add_option('--mpithreads',
                  help = 'Number of threads to start mpi with (default 4)'
                  ,type='int',default=4)
    parser.add_option('--display',
                  help = 'Display image and load catalog on ds9?'
                  ,action='store_true',default=False)
    parser.add_option('--overwrite',
                  help = 'Run in overwrite mode?'
                  ,action='store_true',default=False)

    opt, args = parser.parse_args(argv)

    # Run sextractor on image to get catalog list
    sex = SExtractor()

    # default params
    sex.config['PIXEL_SCALE'] = 0.55
    sex.config['BACK_TYPE'] = "AUTO"
    sex.config['SATUR_LEVEL'] = 60000
    sex.config['DETECT_MINAREA'] = 200
    sex.config['DETECT_THRESH'] = 10.0
    sex.config['VERBOSE_TYPE'] = "QUIET"
    sex.config['PARAMETERS_LIST'] = ["NUMBER",'X_IMAGE', 'Y_IMAGE',
                                     "XWIN_IMAGE", "YWIN_IMAGE",
                                     "FLUX_BEST", "FWHM_IMAGE",
                                     "FLAGS"]

    outname = os.path.basename(opt.image).split('-')[1].replace('.fits','_zern.npy')

    # ok, here we go!
    if ( (not os.path.exists(outname)) or (opt.overwrite) ):

        log.info('Running sextractor')
        sex.run(opt.image, clean=False)

        if opt.display:
            log.info('Displaying image in ds9')

        # Now calculate zernike coef

        cmd = 'mpirun -np %i python ~/Develop/donut/script/donut -i %s -p %s -c %s -o %s'%(opt.mpithreads,
                                                                    opt.image,
                                                                    opt.parameters,
                                                                    sex.config['CATALOG_NAME'],
                                                                    outname)

        log.info('Running donut with %i cores'%opt.mpithreads)
        log.debug(cmd)
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        p.wait()
    else:
        log.warning('File %s exists. Run with --overwrite to force reprocess.'%outname)

    log.info('Mapping hexapod position')

    # zerpar = [os.path.basename(__file__),'-f',outname,'--max', '100', '--niter', '3','-i',opt.image]
    #
    # zernmap.main(zerpar)

    zmap = ZernMap(cfp = 291.36*units.mm,
                   pix2mm = 9*units.micron,
                   center = [9216/2,9232/2]) #[0,0])

    rcat = np.load(outname).T

    cat = np.array([])

    for col in rcat:
        cat = np.append(cat,np.array([col.reshape(-1),]))
    cat = cat.reshape(rcat.shape[0],-1)

    fitmask = cat[0] == 1
    cat = cat[1:]

    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commay = 8
    id_commax = 9

    center = [9216/2,9232/2]

    planeU = zmap.astigmatism(cat[0],cat[1],cat[id_astigx]*sex.config['PIXEL_SCALE'],0)
    planeV = zmap.astigmatism(cat[0],cat[1],cat[id_astigy]*sex.config['PIXEL_SCALE'],1)
    # print cat[id_commax]
    comaX,mask = zmap.comma(cat[0],cat[id_commax]*sex.config['PIXEL_SCALE'])
    comaY,mask = zmap.comma(cat[1],cat[id_commay]*sex.config['PIXEL_SCALE'])
    focus = zmap.map(cat[0],cat[1],cat[id_focus]*sex.config['PIXEL_SCALE'])
    seeing = zmap.map(cat[0],cat[1],cat[id_seeing]*sex.config['PIXEL_SCALE'])

    U = (planeU['U']+planeV['U'])/2.
    V = (planeU['V']+planeV['V'])/2.
    newFitX = np.poly1d(comaX)
    newFitY = np.poly1d(comaY)
    signX = +1
    signY = -1

    print '#'*39

    print '# Offset X: %+6.4f (%+6.4f/%+6.4f) #'%(signX*newFitX(0.)+(planeU['X'].to(units.mm).value+planeV['X'].to(units.mm).value)/2.,signX*newFitX(0.),(planeU['X'].to(units.mm).value+planeV['X'].to(units.mm).value)/2.)
    print '# Offset Y: %+6.4f (%+6.4f/%+6.4f) #'%(signY*newFitY(0.)+(planeU['Y'].to(units.mm).value+planeV['Y'].to(units.mm).value)/2.,signY*newFitY(0.),(planeU['Y'].to(units.mm).value+planeV['Y'].to(units.mm).value)/2.)
    print '# Offset Z: %+6.4f %s #'%(focus[2]/10.,' '*17)
    print '# Offset U: %25s #'%(U.to_string(unit=units.degree, sep=(':', ':', ' ')))
    print '# Offset V: %25s #'%(V.to_string(unit=units.degree, sep=(':', ':', ' ')))


    print '#'*56
    hdr = fits.getheader(opt.image)
    DXHEX = 'HIERARCH T80S TEL FOCU HEX DX'
    DYHEX = 'HIERARCH T80S TEL FOCU HEX DY'
    DZHEX = 'HIERARCH T80S TEL FOCU HEX DZ'
    DUHEX = 'HIERARCH T80S TEL FOCU HEX DU'
    DVHEX = 'HIERARCH T80S TEL FOCU HEX DV'
    print '# Offset X: %+6.4f%+6.4f = %+6.4f (%+6.4f/%+6.4f) #'%(float(hdr[DXHEX]),signX*newFitX(0.)+(planeU['X'].to(units.mm).value+planeV['X'].to(units.mm).value)/2.,
                                                          float(hdr[DXHEX])+signX*newFitX(0.)+(planeU['X'].to(units.mm).value+planeV['X'].to(units.mm).value)/2.,
                                                          float(hdr[DXHEX])+signX*newFitX(0.),float(hdr[DXHEX])+(planeU['X'].to(units.mm).value+planeV['X'].to(units.mm).value)/2.)
    print '# Offset Y: %+6.4f%+6.4f = %+6.4f (%+6.4f/%+6.4f) #'%(float(hdr[DYHEX]),signY*newFitY(0.)+(planeU['Y'].to(units.mm).value+planeV['Y'].to(units.mm).value)/2.,
                                                                 float(hdr[DYHEX])+signY*newFitY(0.)+(planeU['Y'].to(units.mm).value+planeV['Y'].to(units.mm).value)/2.,
                                                                 float(hdr[DYHEX])+signY*newFitY(0.),float(hdr[DYHEX])+(planeU['Y'].to(units.mm).value+planeV['Y'].to(units.mm).value)/2.)
    print '# Offset Z: %+6.4f%+6.4f = %+6.4f %s #'%(float(hdr[DZHEX]),focus[2]/10.,float(hdr[DZHEX])+(focus[2]/10.),' '*17)

    du = Angle(float(hdr[DUHEX])*units.degree)
    dv = Angle(float(hdr[DVHEX])*units.degree)
    corrU = du+U
    corrV = dv+V
    print '# Offset U: %s%s = %s    #'%(du.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                      U.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                      corrU.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True))
    print '# Offset V: %s%s = %s    #'%(dv.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                      V.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True),
                                      corrV.to_string(unit=units.degree, sep=':',precision=2,alwayssign=True,pad=True))
    print '#'*56

    return 0

if __name__ == '__main__':
    main(sys.argv)
