
import numpy as np
from astropy import units
from astropy.coordinates import Angle
from scipy.optimize import leastsq
import logging
#import chimera.core.log
log = logging.getLogger(__name__)

class ZernMapException(Exception):
    pass

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

def mkMask(data,size,nx=3,ny=3):

    grid = np.zeros((nx,ny))

    xgrid = np.arange(size[0])
    ygrid = np.arange(size[1])

    return

class ZernMap:

    @units.quantity_input(cfp=units.meter,pix2mm=units.meter)
    def __init__(self,cfp,pix2mm,center):
        self._CFP = cfp
        self.pix2mm = pix2mm
        self.center = center
        self.maxreject = 100

    @units.quantity_input(val=units.meter)
    def setCFP(self,val):
        self._CFP = val

    def getCFP(self):
        return self._CFP

    def astigmatism(self,x,y,z,Axis=0):

        median = np.median(z)
        std = np.std(z)
        mask = np.abs(z-median) < std*2
        # print median,std
        fitx = (x[mask]-self.center[0])*self.pix2mm.to(units.mm).value
        fity = (y[mask]-self.center[1])*self.pix2mm.to(units.mm).value

        plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)
        ix = Angle(plane[0]*units.rad)/10.
        iy = Angle(plane[1]*units.rad)/10.
        CFP = self.getCFP().to(units.mm)

        zcfp = - CFP * (2. - np.cos(ix.rad) - np.cos(iy.rad))
        xcfp = CFP * np.sin(iy.rad)
        ycfp = CFP * np.sin(ix.rad)

        res = dict([   ('X',0),
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

        return res

    def comma(self,x,y,Axis,xrange=10.,niter=3):

        mask = np.zeros_like(x) == 0

        pos =  (x[mask]-self.center[Axis])*self.pix2mm.to(units.mm).value

        mask = np.bitwise_and(pos > -xrange ,
                              pos <  xrange  )

        zero_reject = len(mask)-len(mask[mask])
        nreject = len(mask)-len(mask[mask])

        pol = None
        for iter in range(niter if niter > 0 else 1):
            x1,y1 = pos[mask],y[mask]

            pol = np.polyfit(x1,y1,1)

            fit = np.poly1d(pol)

            yres = y - fit(pos) # residue
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

            if new_nreject-zero_reject > self.maxreject:
                log.debug('Maximum reject (%i) reached (%i). Breaking.'%(self.maxreject,new_nreject))
                break
            elif new_nreject == nreject:
                log.debug('Rejected 0 data points. Breaking')
                break
            nreject = new_nreject

        if pol is None:
            raise ZernMapException("Could not get comma values.")

        return pol,mask

    def map(self,x,y,z):

        mean = np.mean(z)
        std = np.std(z)
        mask = np.bitwise_and(z < (mean+3*std),
                              z > (mean-3*std))

        if len(mask[mask]) < 100:
            mask = np.zeros_like(z) == 0

        new_x = np.linspace(x.min(),x.max(),101)
        new_y = np.linspace(y.min(),y.max(),101)

        Xi,Yi = np.meshgrid(new_x,new_y)

        fitx = (x[mask]-self.center[0])*self.pix2mm.to(units.mm).value
        fity = (y[mask]-self.center[1])*self.pix2mm.to(units.mm).value

        plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)

        return plane