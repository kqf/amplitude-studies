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

    def getcoef(self, process, odd):
        signature = -1 if process % 2 == 0 else 1
        coef = signature * 1j if odd < 0 else -1
        return coef

    def setup(self, par, i, process = 110):
        self.alpha0, self.alpha, g, self.B, odd = par[i: i + self.npars]
        self.g = g * self.getcoef(process, odd)
        return self.npars

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
        cdllp = log(-1j * s)
        rf2 = self.B + self.alpha * cdllp 
        return self.g * exp(self.alpha0 * cdllp) * (0.5 / rf2) * exp(-0.25 * x ** 2 / rf2) / (8 * pi * s)

    @staticmethod
    def h_amplitude(poles, s, x):
        # print 'len', poles
        return sum(pole.h_(s, x) for i, pole in enumerate(poles))


class TripleExponentAmplitude(Amplitude):
    npars = 9

    def __init__(self):
        super(TripleExponentAmplitude, self).__init__()


    def setup(self, par, i, process):
        delp1, self.alp1p, g, self.betap1, self.betap2, self.betap3, self.cp11, self.cp12, odd = par[i: i + self.npars]
        self.alp1 = par[9] - delp1 if i != 9 else delp1
        self.cp13 = g 
        self.coef = self.getcoef(process, odd)
        return self.npars

    def h_(self, s, x):
        ss = -1j * s
        alphalog = log(ss) * self.alp1p

        f = lambda r: r + alphalog 
        r1, r2, r3 = map(f, (self.betap1, self.betap2, self.betap3))

        fbp1 =(self.cp11 * exp(-0.25 * x **2 / r1) * (0.5 / r1) + 
               self.cp12 * exp(-0.25 * x **2 / r2) * (0.5 / r2) + 
               self.cp13 * exp(-0.25 * x **2 / r3) * (0.5 / r3))

        res = fbp1 * ss ** self.alp1 * self.coef
        return res / (8 * pi * s)



class Model(object):
    model_types = {
        'analytic': (AnalyticAmplitude, [AnalyticAmplitude(), AnalyticAmplitude(), AnalyticAmplitude()]), 
        'numeric': (Amplitude, [Amplitude(), Amplitude(), Amplitude()]), 
        'triples': (AnalyticAmplitude, [TripleExponentAmplitude(), TripleExponentAmplitude(), TripleExponentAmplitude(), AnalyticAmplitude(), AnalyticAmplitude()]),
        'triples_test': (AnalyticAmplitude, [TripleExponentAmplitude()])#, TripleExponentAmplitude(), TripleExponentAmplitude(), AnalyticAmplitude(), AnalyticAmplitude()])
        }
    k, m_p =  0.3893797, 0.9382700

    def __init__(self, mtype = 'analytic', npoles = 3):
        super(Model, self).__init__()
        self.AmplitudeType, self.poles = self.model_types.get(mtype, None)
        self.llambda = 0
        self.observables = {
                             110: self.cross_section, 111: self.cross_section,
                             210: self.rho, 211: self.rho ,
                             310: self.diff_cross_section, 311: self.diff_cross_section
                           }

    def set_parameters(self, code, par):
        i, par = 0, [p for p in par] # ROOT provides buffer, not array
        for pole in self.poles:
            i += pole.setup(par, i, code)
        self.llambda = par[-1]
        
    def __call__(self, s, t, code, par):
        function = self.observables.get(code)
        self.set_parameters(code, par)
        return function(s, t)

    def H_amplitude(self, s, x):
        h = self.AmplitudeType.h_amplitude(self.poles, s, x)
        H = self.unitarize(h)
        return H

    def unitarize(self, h):
        L = complex(0, 2 * self.llambda)
        H = (exp(L * h) - 1) / L
        return H

    def A_amplitude(self, s, t, skipreal = True):
        factor = 8 * pi * s 

        f = lambda b: self.H_amplitude(s, b).imag * b * sf.j0(b * t ** 0.5)
        imag = integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

        if skipreal:
            return imag * factor

        f = lambda b: self.H_amplitude(s, b).real * b * sf.j0(b * t ** 0.5)
        real = integrate.quad(f, 0, 40)[0]  # integral  zero to lower 
        return complex(real, imag) * factor


    def rho(self, s, t = 0):
        A = self.A_amplitude(s, t, False)
        return A.real / A.imag


    def cross_section(self, s, t = 0):
        # flux = sqrt((s - 4 * self.m_p ** 2) * s)
        sigma_tot = (self.k / s) * self.A_amplitude(s, t)
        return sigma_tot


    def diff_cross_section(self, s, t = 0):
        A = self.A_amplitude(s, t, False)
        res = abs(A) ** 2 * self.k / 16./ pi / s  ** 2
        return res


class UMatrixModel(Model):

    def __init__(self):
        super(UMatrixModel, self).__init__()
        pass

    def unitarize(self, h):
        L = complex(0, 2 * self.llambda)
        H = h / ( -0.5 * L * h + complex(1, 0) )
        return H

