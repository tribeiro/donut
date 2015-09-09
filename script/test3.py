
import sys,os
import pylab as py
import numpy as np
from scipy.optimize import leastsq
from astropy import units as u
from astropy.coordinates import Angle

def main(argv):

    _path = '/Volumes/TIAGOSD2/Documents/T80/data/20150831'
    # _path = '/home/tiago/Documents/data/T80S/20150831'
    # files = ['042152_zern.npy','042958_zern.npy']
    # files = ['042958_zern.npy','043232_zern.npy']
    files = ['../20150903/065348_zern.npy','../20150907/063427_zern.npy']
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

    def fitParabola2DOptimize(XYZ,p0=np.array([1.,1.,1.,1.,1.])):
        def residiuals(parameter,f,x,y):
            return [(f[i] - model(parameter,x[i],y[i])) for i in range(len(f))]


        def model(parameter, x, y):
            ax, ay, bx, by, c = parameter
            return ax*x**2 + bx*x + ay*y**2 +by*y+ c

        X = XYZ[:,0]
        Y = XYZ[:,1]
        Z = XYZ[:,2]

        result = leastsq(residiuals, p0, args=(Z,X,Y))[0]

        return result
        normal = result[0:5]
        nn = np.linalg.norm(normal)
        normal = normal / nn
        return normal

    def map():

        median = np.median(z)
        std = np.std(z)
        mask = np.abs(z-median) < std*2
        # print median,std
        fitx = (x[mask]-center[0])*pix2mm
        fity = (y[mask]-center[1])*pix2mm
        plane = fitPlaneOptimize(np.array([fitx,fity,z[mask]]).T)
        ix = Angle(plane[0]*u.rad)/10.
        iy = Angle(plane[1]*u.rad)/10.

        print "Inclination X: %s"%ix.to_string(unit=u.degree, sep=(':', ':', ' '))
        print "Inclination Y: %s"%iy.to_string(unit=u.degree, sep=(':', ':', ' '))
        CFP = 291.36 # Comma Free Point in mm
        zcfp = - CFP * (2. - np.cos(ix.rad) - np.cos(iy.rad))
        xcfp = CFP * np.sin(iy.rad)
        ycfp = CFP * np.sin(ix.rad)
        print 'Corrections:'
        print 'X: ',xcfp
        print 'Y: ',ycfp
        print 'Z: ',zcfp
        newx = np.linspace(-center[0]*pix2mm,center[0]*pix2mm,101)
        newy = np.linspace(-center[1]*pix2mm,center[1]*pix2mm,101)
        XX,YY = np.meshgrid(newx,newy)
        ZZ = plane[0]*XX + plane[1]*YY + plane[2]
        py.pcolor(XX/pix2mm+center[0],YY/pix2mm+center[1],ZZ,vmin=zmin,vmax=zmax)
        # py.colorbar()
        py.scatter(x[mask],y[mask],50,z[mask],marker='o',vmin=zmin,vmax=zmax)
        py.xlim(0,9216)
        py.ylim(0,9232)
        py.colorbar()

        return plane

    def map2(plane1,plane2):

        fitx = (x-center[0])*pix2mm
        fity = (y-center[1])*pix2mm
        parabola = fitParabola2DOptimize(np.array([fitx,fity,z]).T)
        # print parabola
        newx = np.linspace(-center[0]*pix2mm,center[0]*pix2mm,101)
        newy = np.linspace(-center[1]*pix2mm,center[1]*pix2mm,101)
        XX,YY = np.meshgrid(newx,newy)
        ZZ = parabola[0]*XX**2.+parabola[2]*XX+parabola[1]*YY**2+parabola[3]*YY+parabola[4]

        ZZ1 = plane1[0]*XX + plane1[1]*YY + plane1[2]
        ZZ2 = plane2[0]*XX + plane2[1]*YY + plane2[2]
        pXX = XX/pix2mm+center[0]
        pYY = YY/pix2mm+center[1]
        # ZZ = np.sqrt(ZZ1**2+ZZ2**2)
        py.pcolor(pXX,pYY,ZZ,vmin=zmin,vmax=zmax)
        minpix = np.unravel_index(ZZ.argmin(),ZZ.shape)

        cx = -parabola[2]/parabola[0]/2/pix2mm+center[0]
        cy = -parabola[3]/parabola[1]/2/pix2mm+center[1]

        print 'Center @ %fx%f'%(cx,cy)
        print 'Min %.3f @ %fx%f'%(ZZ[minpix[0]][minpix[1]],pXX[minpix[0]][minpix[1]],pYY[minpix[0]][minpix[1]])
        # py.pcolor(XX/pix2mm+center[0],YY/pix2mm+center[1],ZZ,vmin=zmin,vmax=zmax)
        # py.colorbar()
        py.scatter(x,y,50,z,marker='o',vmin=zmin,vmax=zmax)
        py.plot([0,9216],[cy,cy],'k-')
        py.plot([cx,cx],[0,9232],'k-')
        py.xlim(0,9216)
        py.ylim(0,9232)
        py.colorbar()

    zer = np.load(os.path.join(_path,files[0])).T
    zer[2:]/=10
    zmin,zmax = -0.05,0.05


    center = [9216/2,9232/2]
    x = zer[0]
    y = zer[1]
    z = zer[id_astigx]

    py.subplot(231)
    planeU = map()

    x = zer[0]
    y = zer[1]
    z = zer[id_astigy]

    py.subplot(232)
    planeV = map()

    x = zer[0]
    y = zer[1]
    z = np.sqrt(zer[id_astigy]**2.+zer[id_astigx]**2)
    zmin,zmax = 0.0,0.06

    py.subplot(233)
    map2(planeU,planeV)

    zer = np.load(os.path.join(_path,files[1])).T
    #zer[2:]/=10
    zmin,zmax = -0.05,0.05

    x = zer[0]
    y = zer[1]
    z = zer[id_astigx]


    py.subplot(234)
    planeU = map()

    x = zer[0]
    y = zer[1]
    z = zer[id_astigy]

    py.subplot(235)
    planeV = map()

    x = zer[0]
    y = zer[1]
    z = np.sqrt(zer[id_astigy]**2.+zer[id_astigx]**2)
    zmin,zmax = 0.0,0.06

    py.subplot(236)
    map2(planeU,planeV)

    py.show()

    return

if __name__ == '__main__':
    main(sys.argv)