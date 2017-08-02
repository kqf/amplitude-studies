#!/usr/bin/python2.7 

from scipy import integrate
from scipy.special import j0, j1
from math import sqrt
from cmath import exp, log, pi

import numpy
import scipy.special as sf


class Model(object):
    def __init__(self):
        super(Model, self).__init__()
        
    def __call__(self, x, par):
        k =   0.3893797
        m_p = 0.9382700

        s = x[0] ** 2
        flux = sqrt( ( s - 2 * m_p ** 2 ) ** 2 - 4 * m_p ** 4 )
        sigma_tot = (k / flux) * (8 * pi * s ) * self.A_amplitude(s, par)

        print 'sigma total = ' , sigma_tot
        return sigma_tot

    def partial_amplitude(self, s, t, p):
        coef = -1 if p[0] > 0 else 1j
        return p[0] * (-1j*s) ** (p[1] - p[2]*t) * exp(- p[3]*t) * coef

    def a_amplitude(self, s, t,  par):
        return sum( [ self.partial_amplitude(s, t, p) for p in par ] )

    def H_amplitude(self, x, par):
        # Explicit coefitients assignment:
        # TODO: THIS MUST BE REPLACED:
        # TODO: How to pass pjarameters properly ???
        s, parameters = self.rearrange_parameters(par)
        f_real = lambda q: self.a_amplitude(s, q*q, parameters).real * q * sf.j0(q * x)
        h_real = (1. /( 8 * pi * s ) ) * integrate.quad(f_real, 0, numpy.inf)[0]

        f_imag = lambda q: self.a_amplitude(s, q*q, parameters).imag * q * sf.j0(q * x)
        h_imag = (1. /( 8 * pi * s ) ) * integrate.quad(f_imag, 0, numpy.inf)[0]

        h = complex(h_real, h_imag)

        lamda = par[9]
        L = complex(0, 2 * lamda)

        H = h / ( -0.5 * L * h + complex(1, 0) )
        return H

    def A_amplitude(self, s, par):
        parameters = [p for p in par]

        parameters.append(s)

        # For simga_{tot} -- t is unecessary. It is placed here as reminder
        # t = 0
        # parameters.append(t)


        # f = lambda b: H_amplitude(b, parameters).imag * b  
        f = lambda b: self.H_amplitude(b, parameters).imag * b

        return integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

    def rearrange_parameters(self, par):
        return par[-1], [par[0:4], par[4:8], par[8:12]]

        
class AnalyticalModel(Model):

    def __init__(self):
        super(AnalyticalModel, self).__init__()
        pass

    def H_amplitude(self, x, par):
        g_p, alpha_p0, alpha_p, Bp, g_f, alpha_f0, alpha_f, Bf, g_w, alpha_w0, alpha_w, Bw, lamda, s = par
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

        return H