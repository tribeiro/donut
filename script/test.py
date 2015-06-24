
import pylab as py
from donut.don11 import Donut
from astropy.io import fits as pyfits

don = Donut()

don.readpar('/home/tiago/Develop/donut/donut/donut.json')
don.init()

img = pyfits.getdata('/home/tiago/Develop/donut/code/focus.fits')

piximg = don.extract(img.T)

#print piximg.shape
#print don.zres
#x2,immod = don.fit(piximg.T,0)

z = [0.724 , -2.606 , -3.640 ,  1.906 , -0.058 ,  0.243 ,  0.315 , -0.185 , -0.163 ,  0.001  ,-0.125 , -0.095 , -0.061 , -0.020 , -0.029 ,  0.184 , -0.084 , -0.009  , 0.079 , -0.001 , -0.015]
immod = don.getimage(z)

py.subplot(121)
py.imshow(piximg.T,aspect=1,interpolation='nearest',origin='lower')
#
py.subplot(122)
py.imshow(immod,aspect=1,interpolation='nearest',origin='lower')
#
py.show()