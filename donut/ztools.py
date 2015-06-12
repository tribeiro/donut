'''
Translation from IDL to python of donut program developed by Andrei Tokovinin.

 ztools.pro -> ztools.py

 Toolbox of Zernike polynomials
'''

import numpy as np
from scipy.special import gamma

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
            nn,mm=zern_num(jp)
            k_zz = 2.242434
            K = k_zz * (-1)**((n+nn-2.0*m)/2.0) *np.sqrt((n+1.0)*(nn+1.0))

            if (m == mm) and ((j*jp/2.0 != np.ceil(j*jp/2.0)) or ((j/2.0 == np.ceil(j/2.0)) and (jp/2.0 == np.ceil(jp/2.0))) or m == 0):
                c[j-2,jp-2] = K * gamma((n+nn-5.0/3.0)/2.0) / (gamma((n-nn+17.0/3.0)/2.0) * gamma((nn-n+17.0/3.0)/2.0) * gamma((n+nn+23.0/3.0)/2.0))
            else:
                c[j-2,jp-2] = 0.0

    return c


def zern_num(num):
    return 0,0