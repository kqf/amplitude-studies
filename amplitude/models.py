#!/usr/bin/python2.7 

from amplitude import AnalyticAmplitude, NumericAmplitude, TripleExponentAmplitude, NonlinearAmplitude

from cmath import exp, pi
from scipy import integrate
import scipy.special as sf

class Eikonal(object):
    ampl_types = {
                    'analytic': AnalyticAmplitude, 
                    'numeric': NumericAmplitude, 
                    'triples': TripleExponentAmplitude,
                    'nmod':    NonlinearAmplitude
                }

    k, m_p =  0.3893797, 0.9382700

    def __init__(self, mtype = 'analytic'):
        super(Eikonal, self).__init__()
        self.amplitude = self.ampl_types.get(mtype)()
        self.llambda = 0
        self.observables = {
                             110: self.cross_section, 111: self.cross_section,
                             210: self.rho, 211: self.rho ,
                             310: self.diff_cross_section, 311: self.diff_cross_section
                           }

    def set_parameters(self, par):
        self.llambda = self.amplitude.set_parameters(par)

    def __call__(self, s, t, code, par):
        self.set_parameters(par)
        function = self.observables.get(code)
        return function(s, t, code)

    def H_amplitude(self, s, x, code):
        h = self.amplitude.h(s, x, code)
        H = self.unitarize(h)
        return H

    def unitarize(self, h):
        L = complex(0, 2 * self.llambda)
        H = (exp(L * h) - 1) / L
        return H

    def A_amplitude(self, s, t, code, skipreal = True):
        factor = 8 * pi * s 

        f = lambda b: self.H_amplitude(s, b, code).imag * b * sf.j0(b * t ** 0.5)
        imag = integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

        if skipreal:
            return imag * factor

        f = lambda b: self.H_amplitude(s, b, code).real * b * sf.j0(b * t ** 0.5)
        real = integrate.quad(f, 0, 40)[0]  # integral  zero to lower 
        return complex(real, imag) * factor


    def rho(self, s, t = 0, code = 110):
        A = self.A_amplitude(s, t, code, False)
        return A.real / A.imag


    def cross_section(self, s, t = 0, code = 110):
        # flux = sqrt((s - 4 * self.m_p ** 2) * s)
        sigma_tot = (self.k / s) * self.A_amplitude(s, t, code)
        return sigma_tot


    def diff_cross_section(self, s, t = 0, code = 110):
        A = self.A_amplitude(s, t, code, False)
        res = abs(A) ** 2 * self.k / 16./ pi / s  ** 2
        return res


class UMatrix(Eikonal):
    def __init__(self, mtype = 'analytic'):
        super(UMatrix, self).__init__()

    def unitarize(self, h):
        L = complex(0, 2 * self.llambda)
        H = h / ( -0.5 * L * h + complex(1, 0) )
        return H
