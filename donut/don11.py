
import json

'''
Measuring low-order aberrations from defocused images
May 3, 2006
Translation to python Jun, 2015
'''

class Donut():

    def __init__(self):

        # Stores "common block" from IDL routine
        self.donpar = None
        self.ngrid = None
        self.r = None
        self.zgrid = None
        self.pupil = None
        self.Rpix = None
        self.inside = None
        self.asperpix = None
        self.ccdpix = None
        self.npixperpix = None
        self.fovpix = None
        self.sflag = None
        self.sigpix = None
        self.flux = None

    def readpar(self,filename):

        with open(filename) as data_file:
            data = json(data_file)
        