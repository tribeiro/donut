
import os
import pylab as py
from donut.don11 import Donut
from astropy.io import fits as pyfits
import numpy as np

don = Donut()

print os.path.expanduser('~/Develop/WFS_AndreiTokovinin/donut/donut.json')
don.readpar(os.path.expanduser('~/Develop/WFS_AndreiTokovinin/donut/donut.json'))
don.init()

img = pyfits.getdata(os.path.expanduser('~/Develop/WFS_AndreiTokovinin/code/focus.fits'))

piximg = don.extract(img.T)

#print piximg.shape
#print don.zres
x2,immod = don.fit(piximg.T,0)

# z = [0.500   ,-0.315   ,-0.566    ,2.093  , -0.103    ,0.141    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000    ,0.000]
z = [0.724 , -2.606 , -3.640 ,  1.906 , -0.058 ,  0.243 ,  0.315 , -0.185 , -0.163 ,  0.001  ,-0.125 , -0.095 , -0.061 , -0.020 , -0.029 ,  0.184 , -0.084 , -0.009  , 0.079 , -0.001 , -0.015]
#um  0.539   -2.609   -3.326    1.747   -0.085    0.156    0.181   -0.123   -0.125   -0.010    0.089   -0.026   -0.009   -0.018   -0.019    0.039   -0.046   -0.017    0.032   -0.001   -0.011


immod2 = don.getimage(z)

py.subplot(221)
py.imshow(piximg.T,aspect=1,interpolation='nearest',origin='lower')
#
py.subplot(222)
py.imshow(immod2,aspect=1,interpolation='nearest',origin='lower')
#
py.subplot(223)
py.imshow(np.abs(piximg.T-immod),aspect=1,interpolation='nearest',origin='lower')
#
py.subplot(224)
py.imshow(np.abs(piximg.T-immod2),aspect=1,interpolation='nearest',origin='lower')
#
py.show()