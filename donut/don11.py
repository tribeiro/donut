
import numpy as np
import json
import ztools

'''
Measuring low-order aberrations from defocused images
May 3, 2006
Translation to python Jun, 2015
'''

class Donut():

    def __init__(self):
        pass
        # Stores "common block" from IDL routine
        # self.donpar = None
        # self.ngrid = None
        # self.r = None
        # self.zgrid = None
        # self.pupil = None
        # self.Rpix = None
        # self.inside = None
        # self.asperpix = None
        # self.ccdpix = None
        # self.npixperpix = None
        # self.fovpix = None
        # self.sflag = None
        # self.sigpix = None
        # self.flux = None

    def readpar(self,filename):

        with open(filename) as data_file:
            data = json.load(data_file)

        for key in data['donpar']:
            k = key.keys()[0]
            self.__dict__[k.lower()] = key[k]

    def init(self):
        '''
        ; Pre-compute the parameters and save them in the COMMOM block
        :return:
        '''

        ngrid = self.ngrid
        d = self.d
        eps = self.eps
        alambda = self.alambda
        pixel = self.pixel
        self.sflag = 0 # fast seeing blur, set to 1 for slow calc.


        # Find reasonable limits on the grid parameters
        asperpix =  206265.*(1e-6*alambda)/D # maximum fine-pixel size

        ccdpix = float(pixel)
        k = np.floor(np.log10(ccdpix/asperpix)/np.log10(2.)) +1.
        npixperpix = 2**k
        fovpix = 2*ngrid/npixperpix     # CCD field size
        asperpix = ccdpix/npixperpix
        size = 206266.*(1e-6*alambda)/asperpix
        Rpix = ngrid/size*D

        print 'Rebinning factor: ', npixperpix
        print 'Grid pixel: ',asperpix,' arcsec'
        print 'Grid size: ', size,' m'
        print 'CCD pixel:  ',ccdpix,' arcsec'
        print 'CCD field:  ',fovpix*ccdpix,' arcsec'
        print 'CCD format: ', fovpix

        r = np.rool(np.roll(ztools.dist(2*ngrid),
                            ngrid,
                            axis=0),
                    ngrid,
                    axis=1) #distance from grid center, pixs
        inside = np.bitwise_and( r <= Rpix ,
                                 r <= Rpix*eps )
        pupil = np.zeros((2*ngrid,2*ngrid))    # zero array
        pupil[inside] = 1
        n = len(inside)

        x = (np.arange(2*ngrid) - ngrid) # replicate(1.,2*ngrid)
        theta = np.arctan2(x.T,x)
        theta[ngrid][ngrid]=0.
        self.zgrid = np.zeros((n,2))
        self.zgrid[0] = r(inside)/Rpix
        self.zgrid[1] = theta[inside]
        self.inside = inside
        self.fovpix = fovpix
        self.asperpix = asperpix
        self.ccdpix = ccdpix

    def getimage(self,z):
        '''

        :param z: the Zernike vector in microns, starting from Z=2 (tip)
                    z[0] is seeing in arcseconds
        :return:
        '''
        #COMMON imagedata, uampl, filter2, seeing


        fact = 2.*np.pi/self.alambda
        nzer = len(z)
        phase = np.zeros_like(self.zgrid[0]) # empty array for phase
        for j in range(1, nzer):
            phase += fact*z[j]*ztools.zernike_estim(j+1,self.zgrid)

        tmp = np.zeros((self.ngrid*2,self.ngrid*2))
        uampl = np.complex(tmp, tmp)
        self.uampl = uampl
        uampl[self.inside] = np.complex(np.cos(phase), np.sin(phase))

        self.seeing = z[0]

        #---------  compute the image ----------------------
        imh = np.abs(ztools.shift(np.ifft2(ztools.shift(uampl,self.ngrid,self.ngrid)),self.ngrid,self.ngrid))**2.
        if (self.sflag > 0): # exact seeing blur

            filter2 = np.exp(-2.*np.pi**2*(self.seeing/2.35/self.asperpix/2/self.ngrid)**2*self.r**2) # unbinned image
            imh = np.abs(np.fft2(ztools.shift(np.fft2(imh),self.ngrid,self.ngrid)*filter2))
            impix = ztools.rebin(imh,(self.fovpix,self.fovpix)) # rebinning into CCD pixels

        else:
            rr = ztools.shift(ztools.dist(self.fovpix),self.fovpix/2,self.fovpix/2)
            filter2 = np.exp(-2.*np.pi**2*(self.seeing/2.35/self.ccdpix/self.fovpix)**2*rr**2) # binned image
            impix = ztools.rebin(imh,self.fovpix,self.fovpix) # rebinning into CCD pixels
            impix = np.abs(np.fft2(ztools.shift(np.fft2(impix),self.fovpix/2,self.fovpix/2)*filter2)) # Seeing blur


        self.filter2 = filter2
        return impix/np.sum(impix)

    def newimage(self, a, jzer):
        '''
        a is the amplitude change of aberration (micron), Jzer is the Zernike number
        (1=seeing, 4=focus etc.)

        :param a:
        :param jzer:
        :return:
        '''

        #COMMON imagedata, uampl, filter2, seeing

        newampl = self.uampl
        if (jzer > 1): # Change Zernike coefficient
            newphase =  2.*np.pi/self.alambda*a*ztools.zernike_estim(jzer,self.zgrid)
            newampl[self.inside] *= np.complex(np.cos(newphase), np.sin(newphase))
            filter = self.filter2
        else: # new seeing
            newseeing = self.seeing + a
            if (sflag > 0):
                filter = np.exp(-2.*np.pi**2*(newseeing/2.35/self.asperpix/2/self.ngrid)**2*self.r**2) # unbinned image
            else:
                rr = ztools.shift(ztools.dist(self.fovpix),self.fovpix/2,self.fovpix/2)
                filter = np.exp(-2.*np.pi**2*(newseeing/2.35/self.ccdpix/self.fovpix)^2*rr**2) # binned image

        #---------  compute the image ----------------------
        imh = np.abs(ztools.shift(np.ifft2(ztools.shift(newampl,self.ngrid,self.ngrid)),self.ngrid,self.ngrid))**2
        if (sflag > 0): # exact seing blur
            imh = np.abs(np.fft2(ztools.shift(np.fft2(imh),self.ngrid,self.ngrid)*filter))
            impix = ztools.rebin(imh,self.fovpix,self.fovpix) # rebinning into CCD pixels
        else:
            impix = ztools.rebin(imh,self.fovpix,self.fovpix) # rebinning into CCD pixels
            impix = np.abs(np.fft2(ztools.shift(np.fft2(impix),self.fovpix/2,self.fovpix/2)*filter)) # Seeing blur

        return impix/np.sum(impix)

    def getmom(self,impix1): #, impix1, zestim
        '''

        :param impix1:
        :return: Vector of Zernike aberrations in microns
        '''
        n = self.ngrid/self.npixperpix
        xx = (np.arange(2*n)-n)#replicate(1,2*n)
        yy = (np.arange(2*n)-n)#replicate(1,2*n)

        thresh = np.max(impix1)*self.thresh
        impix = impix1
        impix[ impix < thresh ] = 0.

        imh0 = np.sum(impix)

        xc = np.sum(xx*impix)/imh0
        yc = np.sum(yy*impix)/imh0
        mxx = np.sum(impix*(xx-xc)^2)/imh0
        myy = np.sum(impix*(yy-yc)^2)/imh0
        mxy = np.sum(impix*(xx-xc)*(yy-yc))/imh0

        scale = self.npixperpix/(self.ngrid/self.rpix)

        a2 = scale*(xc+0.5)*np.pi*0.5
        a3 = scale*(yc+0.5)*np.pi*0.5
        a4 = scale*np.sqrt((mxx + myy)*0.5)/1.102
        a4 = (a4**2 - (0.5/2.35)**2)
        a4[a4 < 0.] = 0. # subtract 0.5arcsec seeing
        a5 = scale*mxy*(mxx*myy)**(-0.25)/1.45
        a6 = scale*(mxx - myy)*0.5*(mxx*myy)**(-0.25)/1.45
        zestim = [0.,a2,a3,a4,a5,a6]*self.alambda/(2.*np.pi) # estimated Zernike aberrations
        zestim[0] = 0.5

        return zestim
