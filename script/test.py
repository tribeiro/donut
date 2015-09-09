
import os
import numpy as np
import pylab as py
from donut.don11 import Donut
from astropy.io import fits as pyfits
import logging

logging.basicConfig(format='%(asctime)s[%(levelname)s]-%(name)s-(%(filename)s:%(lineno)d):: %(message)s',
                    level=logging.DEBUG)

don = Donut()

# don.readpar(os.path.expanduser('~/Develop/donut/donut/donut.json'))
don.readpar('/Volumes/TIAGOSD2/Documents/T80/data/hexapod/donut.json')
don.init()

# img = pyfits.getdata(os.path.expanduser('~/Develop/donut/code/focus.fits'))
img = pyfits.getdata('/Volumes/TIAGOSD2/Documents/T80/data/hexapod/foo_0029.fits')

piximg = don.extract(img.T)

#print piximg.shape
#print don.zres

#impix_idl = np.loadtxt(os.path.expanduser("~/Develop/donut/code/impix.dat"))
#impix_idl = impix_idl.reshape(16,16)

#immod_idl = np.loadtxt(os.path.expanduser("~/Develop/donut/code/immod.dat")).reshape(16,16)

# x2,immod = don.fit(impix_idl,0)
# x2,immod,mask = don.fit(piximg.T)

z = [0.724 , -2.606 , -3.640 ,  1.906 , -0.058 ,  0.243 ,  0.315 , -0.185 , -0.163 ,  0.001  ,-0.125 , -0.095 , -0.061 , -0.020 , -0.029 ,  0.184 , -0.084 , -0.009  , 0.079 , -0.001 , -0.015]

#um    0.462   -2.647   -3.264    1.753   -0.085    0.150    0.158   -0.107   -0.114   -0.010    0.105   -0.023   -0.007   -0.017   -0.017    0.032   -0.037   -0.013    0.029   -0.000   -0.008
#um    0.569   -2.599   -3.349    1.750   -0.085    0.158    0.190   -0.128   -0.127   -0.008    0.084   -0.028   -0.010   -0.018   -0.020    0.041   -0.049   -0.018    0.034   -0.001   -0.012

#z = [0.724 , -2.606 , -3.640 ,  1.906 , -0.058 ,  0.243 ,  0.315 , -0.185 , -0.163 ,  0.001  ,-0.125 , -0.095 , -0.061 , -0.020 , -0.029 ,  0.184 , -0.084 , -0.009  , 0.079 , -0.001 , -0.015]
# z = np.zeros(len(z))
# z[0] = 0.8
# #z[3] = 1.906
z = np.zeros(11)
z[0] = 1.7
z[3] = 0.
immod = don.getimage(z)

hdu = pyfits.PrimaryHDU(data=immod)
hdu.writeto('/tmp/foo.fits')

py.figure(1)
#
py.subplot(221)
py.imshow(piximg.T,aspect=1,interpolation='nearest',origin='lower')
#
# py.subplot(223)
# py.imshow(impix_idl,aspect=1,interpolation='nearest',origin='lower')
#
py.subplot(222)
py.imshow(immod,aspect=1,interpolation='nearest',origin='lower')
#py.imshow(piximg.T/impix_idl,aspect=1,interpolation='nearest',origin='lower')
#
# py.subplot(224)
# py.imshow(immod_idl,aspect=1,interpolation='nearest',origin='lower')

# #
# py.figure(2)
# #
# py.subplot(221)
# mask = immod < np.mean(immod) *0.25
# comp = (immod-immod_idl)/immod_idl*100.
# comp[mask] = 0.
# # comp[np.bitwise_not(mask)] = 1.
# imax = py.imshow(comp,aspect=1,interpolation='nearest',origin='lower')
# cbar = py.colorbar(imax)
# #
# py.subplot(222)
# comp = (immod-piximg.T)/piximg.T*100.
# #comp[np.bitwise_not(mask)] = 0.
# imax = py.imshow(comp,aspect=1,interpolation='nearest',origin='lower')
# cbar = py.colorbar(imax)
# #
# print np.std(comp)
# py.subplot(223)
# comp = (immod_idl-piximg.T)/piximg.T*100.
# #comp[np.bitwise_not(mask)] = 0.
# imax = py.imshow(comp,aspect=1,interpolation='nearest',origin='lower')
# cbar = py.colorbar(imax)
# #
# print np.std(comp)
py.show()