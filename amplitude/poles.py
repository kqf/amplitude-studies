from scipy import integrate
from scipy.special import j0, j1
from math import sqrt
from cmath import exp, log, pi

import numpy
import scipy.special as sf

class PoleNumeric(object):
    npars = 5
    def __init__(self):
        super(PoleNumeric, self).__init__()

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
        f_real = lambda q: PoleNumeric.a_amplitude(poles, s, q*q).real * q * sf.j0(q * x)
        h_real = (1. / (8 * pi * s) ) * integrate.quad(f_real, 0, numpy.inf)[0]

        f_imag = lambda q: PoleNumeric.a_amplitude(poles, s, q*q).imag * q * sf.j0(q * x)
        h_imag = (1. / (8 * pi * s) ) * integrate.quad(f_imag, 0, numpy.inf)[0]
        return complex(h_real, h_imag)


class Pole(PoleNumeric):
    def __init__(self):
        super(Pole, self).__init__()

    def h_(self, s, x):
        cdllp = log(-1j * s)
        rf2 = self.B + self.alpha * cdllp 
        return self.g * exp(self.alpha0 * cdllp) * (0.5 / rf2) * exp(-0.25 * x ** 2 / rf2) / (8 * pi * s)

    @staticmethod
    def h_amplitude(poles, s, x):
        return sum(poleNumeric.h_(s, x) for poleNumeric in poles)


class TripleExponentPole(PoleNumeric):
    npars = 9

    def __init__(self):
        super(TripleExponentPole, self).__init__()

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

