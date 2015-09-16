#!/usr/env python

'''
Pipeline for determining hexapod corrections.
'''

import sys,os
from chimera.util.sextractor import SExtractor
import logging
import subprocess
import zernmap

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

    # ok, here we go!

    log.info('Running sextractor')
    sex.run(opt.image, clean=False)

    if opt.display:
        log.info('Displaying image in ds9')

    # Now calculate zernike coef

    outname = os.path.basename(opt.image).split('-')[1].replace('.fits','_zern.npy')

    cmd = 'mpirun -np %i python ~/Develop/donut/script/donut -i %s -p %s -c %s -o %s'%(opt.mpithreads,
                                                                opt.image,
                                                                opt.parameters,
                                                                sex.config['CATALOG_NAME'],
                                                                outname)

    if not os.path.exists(outname) and not opt.overwrite:
        log.info('Running donut with %i cores'%opt.mpithreads)
        log.debug(cmd)
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        p.wait()
    else:
        log.warning('File %s exists. Run with --overwrite to force reprocess.'%outname)

    log.info('Mapping hexapod position')

    zerpar = [os.path.basename(__file__),'-f',outname,'--max', '100', '--niter', '3','-i',opt.image]

    zernmap.main(zerpar)

    return 0

if __name__ == '__main__':
    main(sys.argv)
