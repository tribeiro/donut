'''
Translation from IDL to python of donut program developed by Andrei Tokovinin.

 ztools.pro -> ztools.py

 Toolbox of Zernike polynomials
'''

import numpy as np
from scipy.special import gamma,jv
from scipy.linalg import svd

class ZToolsException(Exception):
    pass

def cova_zern1(jmax):
    '''
    Compute the Noll covariance matrix of Zernike polynomials

    :param jmax:
    :return: noll covariance matrix
    '''

    c = np.zeros((jmax-1,jmax-1))

    for j in range(2,jmax):
        n,m = zern_num(j)
        for jp in range(2,jmax):
            nn,mm= zern_num(jp)
            k_zz = 2.242434
            K = k_zz * (-1)**((n+nn-2.0*m)/2.0) *np.sqrt((n+1.0)*(nn+1.0))

            if (m == mm) and ((j*jp/2.0 != np.ceil(j*jp/2.0)) or ((j/2.0 == np.ceil(j/2.0)) and (jp/2.0 == np.ceil(jp/2.0))) or m == 0):
                c[j-2,jp-2] = K * gamma((n+nn-5.0/3.0)/2.0) / (gamma((n-nn+17.0/3.0)/2.0) * gamma((nn-n+17.0/3.0)/2.0) * gamma((n+nn+23.0/3.0)/2.0))
            else:
                c[j-2,jp-2] = 0.0

    return c


def zern_num(num):
    '''
    Description
    ZERN_NUM computes the azimuthal degree M and the radial order N
    for the sequential Zernike polynomial NUM
    INFO  =  result display
    :param num:
    :return: N,M,INFO
    '''

    if num < 1:
        raise ZToolsException('NUM must be an integer greater or equal 1')

    j = float(num)
    n = np.sqrt(8.*j-7.)-1./2.

    if n%2. == 0.:
        # even n
        m = 2.*((j-(n*(n+1.))/2.)/2.)
    else:
        # odd n
        m = 1.+2.*((j-1.-(n*(n+1.))/2.)/2.)

    return n,m

def zernike_estim(mode, grid):
    '''

    :param mode:
    :param grid: set of points (polar co.) on which the Zernike must be evaluated
    :return: vector of Zernike values on the grid
    '''
    n,m = zern_num(mode)
    p = (mode%2)

    R=0.

    for J in range(int((n-m)/2)):
        S= J
        R= R+(-1.)**J*Fact(n-J)/(Fact(S)*Fact((n+m)/2-J)*Fact((n-m)/2-J))*grid[0]**(n-2*J)

    if (m == 0):
        return np.sqrt(n+1.0)*R
    elif (p == 0):
        return np.sqrt(n+1.0)*np.sqrt(2.0)*np.cos(m*grid[1])*R
    else:
        return np.sqrt(n+1.0)*np.sqrt(2.0)*np.sin(m*grid[1])*R

def svd_invert(matrix,threshold):
    '''
    :param matrix:
    :param threshold:
    :return:SCD-inverted matrix
    '''
    u,ws,v = svd(matrix)

    ww = np.max(ws)
    n = len(ws)
    invw = np.identity(n)
    ncount = 0

    for i in range(n):
        if ws[i] < ww*threshold:
            invw[i][i]= 0.
            ncount+=1
        else:
            invw[i][i] = 1./ws[i]

    print '%i singular values rejected in inversion'%ncount

    inv_matrix = np.dot( np.dot(v, invw), np.transpose(u) )

    return inv_matrix

def zern_derivx(j):
    '''
    ; This function calculates the x and y derivative coefficients needed to compute
    ; the derivative of the jth Zernike polynomial.
    ; (d/dx)Zj=SUM_j' gammax_j' Z_j'
    ; (d/dy)Zj=SUM_j' gammay_j' Z_j'
    ; gammax and gammay is the output vector gamma=[2,j]
    ;
    ; Date : 9 December 1999
    ; Written by Elise Viard, eviard@eso.org
    ; Translated to python by Tiago Ribeiro tribeiro@ufs.br

    :param j:
    :return:
    '''

    n,m = zern_num(j)
    gam = np.zeros(j)

    for j2 in range(1,j+1):
        n2,m2 = zern_num(j2)
        if ( (m-m2)**2 == 1):
            if (m != 0) and (m2 != 0):
                if ((j%2 == 0) and (j2%2 == 0)) or ((j%2 != 0) and (j2%2 != 0)):
                    gam[j2-1] = np.sqrt( (n+1)*(n2+1))
                else:
                    gam[j2-1] = 0.
            elif ( (m==0) and (j2%2 == 0)):
                gam[j2-1] = np.sqrt(2.*(n+1)*(n2+1))
            elif ( (m2 ==0) and (j%2 == 0) ):
                gam[j2-1] = np.sqrt(2.*(n+1)*(n2+1))
            else:
                gam[j2-1] = 0
        else:
            gam[j2-1] = 0
    return gam

def zern_derivy(j):
    n,m = zern_num(j)
    gam = np.zeros(j)

    for j2 in range(1,j+1):
        n2,m2 = zern_num(j2)
        if ( (m-m2)**2 == 1):
            if (m != 0) and (m2 != 0):
                if ((j%2 == 0) and (j2%2 == 0)) or ((j%2 != 0) and (j2%2 != 0)):
                    sig = 1.
                    if m2 == m+1 and (j%2 != 0):
                        sig = -1.
                    elif m2 == m-1 and (j%2 == 0):
                        sig = -1.
                    gam[j2-1] = sig*np.sqrt( (n+1)*(n2+1))
                else:
                    gam[j2-1] = 0.
            elif ((m==0) and (j2%2 == 0)):
                gam[j2-1] = np.sqrt(2.*(n+1)*(n2+1))
            elif ((m2 ==0) and (j%2 == 0)):
                gam[j2-1] = np.sqrt(2.*(n+1)*(n2+1))
            else:
                gam[j2-1] = 0
        else:
            gam[j2-1] = 0

    return gam

def zern_deriv(j):
    gam = np.zeros((2,j))
    gam[0] = zern_derivx(j)
    gam[1] = zern_derivy(j)
    return gam

def getftzer(Jzer,ngrid=128,Rpix=100):
    '''
    ; Compute the Fourier Transform of Zernike mode

    ; ngrid = 128 ; grid half-size, pixels
    ; Rpix = 100 ; pupil radius in pixels

    :param Jzer:
    :return:
    '''

    x = np.arange(-ngrid,ngrid)
    y = np.arange(-ngrid,ngrid)
    theta = np.arctan2(x,y)

    n,m = zern_num(Jzer)
    f = np.roll(np.roll(dist(2*ngrid),
                        ngrid,
                        axis=0),
                ngrid,
                axis=1)/(2*ngrid)*Rpix
    f[ngrid][ngrid] = 1e-3

    ftmod = np.sqrt(n+1.0)*jv(n+1,2*np.pi*f)/(np.pi*f)

    if m == 0:
        zz = ftmod*np.complex(0,1.)**(n/2.)
    else:
        if (Jzer%2 == 0):
            fact=np.sqrt(2.)*np.cos(m*theta)
        else:
            fact=np.sqrt(2.)*np.sin(m*theta)
        zz = ftmod*fact*(-1)**((n-m/2.))*np.complex(0,1.)**m

    return zz

def Fact(n):
    return np.math.factorial(n)

def dist(size):
   """
   title::
      dist

   description::
      This method will create a rectangular array (numpy.ndarray) in which
      each element represents the distance from the upper left position -or-
      in which each element is proportional to its frequency.  This method
      is modeled directly after the IDL DIST function.

   attributes::
      size
         a scalar or 2-element tuple/list representing the square dimension
         or the (rows, columns) of the array to be producedi, respectively.

   author::
      Carl Salvaggio

   copyright::
      Copyright (C) 2014, Rochester Institute of Technology

   license::
      GPL

   version::
      1.0.0

   disclaimer::
      This source code is provided "as is" and without warranties as to
      performance or merchantability. The author and/or distributors of
      this source code may have made statements about this source code.
      Any such statements do not constitute warranties and shall not be
      relied on by the user in deciding whether to use this source code.

      This source code is provided without any express or implied warranties
      whatsoever. Because of the diversity of conditions and hardware under
      which this source code may be used, no warranty of fitness for a
      particular purpose is offered. The user is advised to test the source
      code thoroughly before relying on it. The user must assume the entire
      risk of using the source code.
   """

   (rows, columns) = size if isinstance(size, (list, tuple)) else (size, size)

   x = numpy.arange(columns, dtype=numpy.float32)
   x = numpy.where(x < (columns-x), x**2, (columns-x)**2)
   a = numpy.zeros((rows, columns), dtype=numpy.float32)
   for i in range(rows/2+1):
      y = numpy.sqrt(x + i**2)
      a[i,:] = y
      if i != 0:
         a[rows-i,:] = y

   return a


def shift(matrix,s1,s2):
    return np.roll(np.roll(matrix,
                           s1,
                           axis=1),
                   s2,
                   axis=2)

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)