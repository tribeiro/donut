__author__ = 'tiago'

import numpy as np
from donut.don11 import Donut
from tempfile import NamedTemporaryFile
from astropy.io import fits
import json
import logging

log = logging.getLogger(__name__)

class runDonut():
    def __init__(self,parcatalog,image,jobgrid):

        self.parcatalog = parcatalog
        self.donut = Donut()
        self.tmpCat = NamedTemporaryFile(delete=False)
        with open(self.tmpCat.name,'w') as fp:
            json.dump(parcatalog,fp)
        self.img = fits.getdata(image)
        self.jobgrid = jobgrid
        self.fitted = np.zeros_like(jobgrid) == 0

    def donutfit(self):

        log.info('Working on %i jobs (%i:%i)'%(len(self.jobgrid),
                                                    self.jobgrid[0],
                                                    self.jobgrid[-1]))

        self.donut.readpar(self.tmpCat.name,0)
        self.donut.init()
        log.info('Starting ')
        # return
        piximg = self.donut.extract(self.img.T)

        for job in self.jobgrid:
            # don.readpar(tmpCat.name,job)
            # don.init()
            # return
            # img = fits.getdata(opt.image)
            # piximg = don.extract(img.T)

            # for i in range(int(3e5)):
            #     b = np.sum(np.zeros((100,100)))
            # continue

            try:
                if piximg.shape == (self.donut.ngrid/2,self.donut.ngrid/2):
                    # log.info('Fit Start')
                    x2,immod,zer = self.donut.fit(piximg.T)
                    # log.info('Fit End')
                    # zres[job][0] = don.xc
                    # zres[job][1] = don.yc
                    # zres[job][2:] = zer
                else:
                    log.warning('Source %i to close to the border. Skipping...'%(job))
                    fitted[job] = False
            except AttributeError,e:
                log.exception(e)
                fitted[job] = False
                pass
            except Exception,e:
                raise
        log.info('Done')

        return
