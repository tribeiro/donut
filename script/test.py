
from donut.don11 import Donut
from astropy.io import fits as pyfits

don = Donut()

don.readpar('/Users/tiago/Develop/WFS_AndreiTokovinin/donut/donut.json')
don.init()

img = pyfits.getdata('/Users/tiago/Develop/WFS_AndreiTokovinin/code/focus.fits')

piximg = don.extract(img)