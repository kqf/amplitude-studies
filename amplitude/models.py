#!/usr/bin/python2.7 

from scipy import integrate
from scipy.special import j0, j1
from math import sqrt
from cmath import exp, log, pi

import numpy
import scipy.special as sf


class Amplitude(object):
    npars = 5
    def __init__(self):
        super(Amplitude, self).__init__()

    def setup(self, par, i):
        par[self.npars * i: self.npars * (i + 1)] 
        g, self.alpha0, self.alpha, self.B, odd = par[self.npars * i: self.npars * (i + 1)]
        self.odd = odd < 0
        coef = 1j if self.odd else -1
        self.g = g * coef

    def partial_amplitude(self, s, t):
        return self.g * (-1j*s) ** (self.alpha0 - self.alpha * t) * exp(-self.B * t)

    @staticmethod
    def a_amplitude(poles, s, t):
        return sum(p.partial_amplitude(s, t) for p in poles)

    @staticmethod
    def h_amplitude(poles, s, x):
        f_real = lambda q: Amplitude.a_amplitude(poles, s, q*q).real * q * sf.j0(q * x)
        h_real = (1. / (8 * pi * s) ) * integrate.quad(f_real, 0, numpy.inf)[0]

        f_imag = lambda q: Amplitude.a_amplitude(poles, s, q*q).imag * q * sf.j0(q * x)
        h_imag = (1. / (8 * pi * s) ) * integrate.quad(f_imag, 0, numpy.inf)[0]
        return complex(h_real, h_imag)


class AnalyticAmplitude(Amplitude):
    def __init__(self):
        super(AnalyticAmplitude, self).__init__()

    def h_(self, s, x):
        S = -1j * s
        logS_a = log(S ** self.alpha)
        res = ( (self.g / (16. * pi * s))
              * ( ( S ** self.alpha0 ) / ( logS_a + self.B) )
              * exp( -(x / 2.) ** 2  / ( logS_a + self.B) )
             )
        return res

    @staticmethod
    def h_amplitude(poles, s, x):
        return sum(pole.h_(s, x) for pole in poles)

  
class Model(object):
    def __init__(self, mtype = 'analytic', npoles = 3):
        super(Model, self).__init__()
        self.Amplitude = {'analytic': AnalyticAmplitude, 'numeric': Amplitude}.get(mtype, None)
        self.poles = [self.Amplitude() for i in range(npoles)]
        self.llambda = 0

    def set_parameters(self, code, par):
        parameters = [p for p in par]
        for i, pole in enumerate(self.poles):
            pole.setup(parameters, i)
        self.llambda = parameters[-1]
        
    def __call__(self, x, code, par):
        k, m_p, s =  0.3893797, 0.9382700, x ** 2
        self.set_parameters(code, par)

        flux = sqrt( ( s - 2 * m_p ** 2 ) ** 2 - 4 * m_p ** 4 )
        sigma_tot = (k / flux) * (8 * pi * s ) * self.A_amplitude(s)
        return sigma_tot

    def H_amplitude(self, s, x):
        h = self.Amplitude.h_amplitude(self.poles, s, x)
        H = self.unitarize(h)
        return H

    def unitarize(self, h):
        L = complex(0, 2 * self.llambda)
        H = h / ( -0.5 * L * h + complex(1, 0) )
        return H

    def A_amplitude(self, s):
        f = lambda b: self.H_amplitude(s, b).imag * b
        return integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

        
class EikonalModel(Model):

    def __init__(self):
        super(EikonalModel, self).__init__()
        pass

    def unitarize(self, x):
        L = complex(0, 2 * self.llambda)
        H = (exp(L * h) - 1) / L
        return H