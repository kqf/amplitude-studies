#!/usr/bin/python2.7 

from scipy import integrate
from scipy.special import j0, j1
from math import sqrt
from cmath import exp, log, pi

import numpy


N_PARAMETERS_PER_POLE = 4
N_POLES = 3


def fit_function(x, par):
    k =   0.3893797
    m_p = 0.9382700

    s = x[0] ** 2
    flux = sqrt( ( s - 2 * m_p ** 2 ) ** 2 - 4 * m_p ** 4 )
    sigma_tot = (k / flux) * (8 * pi * s ) * A_amplitude(s, par)
    print sigma_tot

    return sigma_tot

def partial_amplitude(s, t, p):
    """ Naiv but cruicial design. Values of parameters
        0  - g
        1  - alpha
        2  - alpha_prime
        3  - B 
    """
    return p[0] * (-1j*s) ** (p[1] - p[2]*t) * exp(- p[3]*t)

def a_amplitude(s, t,  par):
    """ Total gauge amplitude"""
    return sum( [ partial_amplitude(s, t, p) for p in par ] )

def H_amplitude(x, par):
    # This case needs explicit coefitient assignment:
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

def H_amplitude_numerical(x, par):
    """Numerical approach"""
    # Explicit coefitients assignment:
    # TODO: THIS MUST BE REPLACED:
    # TODO: How to pass pjarameters properly ???
    s = par[13]
    # g_p, g_f, g_w, alpha_p, alpha_p0, alpha_f, alpha_f0, alpha_w, alpha_w0, lamda, Bp, Bf, Bw, s, pt = par

    parameters = rearrange_parameters(par)


    f_real = lambda x: a_amplitude(s, x, parameters).real
    h_real = (1. /( 8 * pi * s ) ) * integrate.quad(f_real, 0, numpy.inf)[0]

    f_imag = lambda x: a_amplitude(s, x, parameters).imag
    h_imag = (1. /( 8 * pi * s ) ) * integrate.quad(f_imag, 0, numpy.inf)[0]

    h = complex(h_real, h_imag)

    lamda = par[9]
    L = complex(0, 2 * lamda)

    H = h / ( -0.5 * L * h + complex(1, 0) )
    return H

def A_amplitude(s, par):
    parameters = [ p for p in par ]

    parameters.append(s)

    # For simga_{tot} -- t is unecessary. It is placed here as reminder
    t = 0
    parameters.append(t)


    # f = lambda x: H_amplitude(x, parameters)
    f = lambda x: H_amplitude_numerical(x, parameters).imag

    return integrate.quad(f, 0, numpy.inf)[0]  # integral  zero to lower 

def rearrange_parameters(par):
    """
          Simple convertation from one input parameters ordering convention to another.
          Output parameters presented in for of 2d list:
              [ 
                [g1, a1, ap1, B1],
                ... 

                [gn, an, apn, Bn]]
              ]
    """
    # g_p, g_f, g_w, alpha_p, alpha_p0, alpha_f, alpha_f0, alpha_w, alpha_w0, lamda, Bp, Bf, Bw, s, pt = par
    pomeron = [ par[0], par[4], par[3], par[10] ]
    f_regge = [ par[1], par[6], par[5], par[11] ]
    w_regge = [ par[2], par[8], par[7], par[12] ]
    return [pomeron, f_regge, w_regge]
