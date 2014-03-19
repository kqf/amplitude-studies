#!/usr/bin/python2.7 

from scipy import integrate
from scipy.special import j0, j1
from math import sqrt
from cmath import exp, log, pi
import numpy

def fit_function(x, par):
    k =   0.3893797
    m_p = 0.9382700

    s = x[0] ** 2
    flux = sqrt( ( s - 2 * m_p ** 2 ) ** 2 - 4 * m_p ** 4 )
    sigma_tot = (k / flux) * (8 * pi * s ) * A_amplitude(s, par)

    return sigma_tot

def H_amplitude(x, par):
    g_p, g_f, g_w, alpha_p, alpha_p0, alpha_f, alpha_f0, alpha_w, alpha_w0, lamda, Bp, Bf, Bw, s, pt = par
    alpha  = [alpha_p, alpha_f, alpha_w]
    alpha0 = [alpha_p0, alpha_f0, alpha_w0]
    S = complex(0, -1.*s)
    g1 = complex(-g_p, 0)
    g2 = complex(-g_f, 0)
    g3 = complex(0, -g_w)
    g = [g1, g2, g3]
    B = [Bp, Bf, Bw]
    L = complex(0, 2 * lamda)

    h = complex(0, 0)
    for i in range(3):
        logS_a = log( S ** alpha[i] )
        h += ( (g[i] / (16. * pi * s))
              * ( ( S ** alpha0[i] ) / ( logS_a + B[i] ) )
              * exp( -(x / 2.) ** 2  / ( logS_a + B[i] ) )
             )

    H = h / ( -0.5 * L * h + complex(1, 0) )
    return (H * x).imag


def A_amplitude(s, par):
    parameters = [ p for p in par ]

    parameters.append(s)
    t = 0
    parameters.append(t)

    f = lambda x: H_amplitude(x, parameters)
    return integrate.quad(f, 0, numpy.inf)[0]  # integral  zero to lower 


