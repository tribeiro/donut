
import sys,os
import pylab as py
import numpy as np
from scipy.optimize import leastsq

def main(argv):

    _path = '/Volumes/TIAGOSD2/Documents/T80/data/20150831'
    # files = ['042152_zern.npy','042958_zern.npy']
    # files = ['042958_zern.npy','043232_zern.npy']
    files = ['042152_zern.npy','043616_zern.npy']
    pix2mm = 0.01 # pixel size in um
    id_seeing = 2
    id_focus = 5
    id_astigx = 6
    id_astigy = 7
    id_commay = 8
    id_commax = 9

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
        normal = result[0:3]
        nn = np.linalg.norm(normal)
        normal = normal / nn
        return normal

    def map():

        median = np.median(z)
        std = np.std(z)
        mask = np.abs(z-median) < std*2
        # print median,std
        fitx = x[mask]/np.max(x[mask])-0.5
        fity = y[mask]/np.max(y[mask])-0.5
        plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)
        print plane
        newx = np.linspace(-0.5,0.5,101)
        newy = np.linspace(-0.5,0.5,101)
        XX,YY = np.meshgrid(newx,newy)
        ZZ = plane[0]*XX + plane[1]*YY + plane[2]
        py.pcolor((XX+0.5)*np.max(x[mask]),(YY+0.5)*np.max(y[mask]),ZZ,vmin=zmin,vmax=zmax)
        py.colorbar()
        py.scatter(x[mask],y[mask],50,z[mask],marker='o',vmin=zmin,vmax=zmax)
        py.xlim(0,9216)
        py.ylim(0,9232)
        py.colorbar()

    zer = np.load(os.path.join(_path,files[0])).T
    zmin,zmax = -0.5,0.5


    x = zer[0]
    y = zer[1]
    z = zer[id_astigx]

    py.subplot(231)
    map()

    x = zer[0]
    y = zer[1]
    z = zer[id_astigy]

    py.subplot(232)
    map()

    x = zer[0]
    y = zer[1]
    z = np.sqrt(zer[id_astigy]**2.+zer[id_astigx]**2)
    zmin,zmax = 0.0,0.6

    py.subplot(233)
    map()

    zer = np.load(os.path.join(_path,files[1])).T
    zmin,zmax = -0.5,0.5

    x = zer[0]
    y = zer[1]
    z = zer[id_astigx]

    py.subplot(234)
    map()

    x = zer[0]
    y = zer[1]
    z = zer[id_astigy]

    py.subplot(235)
    map()

    x = zer[0]
    y = zer[1]
    z = np.sqrt(zer[id_astigy]**2.+zer[id_astigx]**2)
    zmin,zmax = 0.0,0.6

    py.subplot(236)
    map()

    py.show()

    return

if __name__ == '__main__':
    main(sys.argv)