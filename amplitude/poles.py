from scipy import integrate
from cmath import exp, log, pi, sqrt

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

    def setup(self, par, i):
        self.alpha0, self.alpha, g, self.B, odd = par[i: i + self.npars]
        self.coef = lambda process: g * self.getcoef(process, odd)
        return self.npars

    def partial_amplitude(self, s, t, process):
        g = (-1j * s) ** (self.alpha0 - self.alpha * t)
        return self.coef(process) * g * exp(-self.B * t)

    @staticmethod
    def a_amplitude(poles, s, t, process):
        return sum(p.partial_amplitude(s, t, process) for p in poles)

    @classmethod
    def h_amplitude(klass, poles, s, x, process):
        def f_real(q):
            return klass.a_amplitude(
                poles, s, q * q, process).real * q * sf.j0(q * x)
        h_real = integrate.quad(f_real, 0, numpy.inf)[0]

        def f_imag(q):
            return klass.a_amplitude(
                poles, s, q * q, process).imag * q * sf.j0(q * x)
        h_imag = integrate.quad(f_imag, 0, numpy.inf)[0]
        return complex(h_real, h_imag) / (8 * pi * s)


class Pole(PoleNumeric):
    def __init__(self):
        super(Pole, self).__init__()

    def h_(self, s, x, process):
        ss = -1j * s
        rf2 = self.B + self.alpha * log(ss)
        hh = ss ** self.alpha0 * (0.5 / rf2) * exp(-0.25 * x ** 2 / rf2)
        return self.coef(process) * hh / (8 * pi * s)

    @classmethod
    def h_amplitude(klass, poles, s, x, process):
        return sum(poleNumeric.h_(s, x, process) for poleNumeric in poles)


class TripleExponentPole(Pole):
    npars = 9

    def __init__(self):
        super(TripleExponentPole, self).__init__()

    def setup(self, par, i):
        delp1, self.alphap, betap1, betap2, betap3, \
            cp11, cp12, cp13, odd = par[i: i + self.npars]

        self.cp = cp11, cp12, cp13
        self.beta = betap1, betap2, betap3
        # TODO: Since this was the bug, try to fix this in next datasets
        #
        self.cp, self.beta = self.beta, self.cp

        self.alpha0 = par[9] - delp1 if i != 9 else delp1
        self.coef = lambda process: self.getcoef(process, odd)
        return self.npars

    def cpsum(self, x, rr):
        return sum(cp * (0.5 / r) * exp(-0.25 * x ** 2 / r)
                   for cp, r in zip(self.cp, rr))

    def h_(self, s, x, process):
        ss = -1j * s
        rr = map(lambda r: r + log(ss) * self.alphap, self.beta)
        fbp1 = self.cpsum(x, rr)

        res = fbp1 * ss ** self.alpha0 * self.coef(process)
        return res / (8 * pi * s)


class NonlinearPole(PoleNumeric):
    npars = 1

    def __init__(self):
        super(NonlinearPole, self).__init__()

    def setup(self, par, i):
        # TODO: Check parameter order
        self.be1, self.be2, self.g1, self.g2, self.g3, \
            self.x01, self.x02, self.amu, self.z, odd = par[i: i + self.npars]

        self.coef = lambda process: self.getcoef(process, odd)
        return self.npars

    def partial_amplitude(self, s, t, process):
        first = (self.g1 * exp(-self.be1 * self.z) +
                 self.g2 * exp(-self.be2 * (
                     sqrt(self.z**2 + self.x02**2) - self.x02)))
        second = self.g3 * (1. + 1. / (1. + self.z / self.x03) ** self.amu)**2
        return first + second


class FirstPomeron(PoleNumeric):
    npars = 10

    def __init__(self):
        super(PoleNumeric, self).__init__()

    def setup(self, par, i):
        delp1, self.alphap, self.t0p, \
            self.g1, self.g2, self.g3, \
            self.beta1, self.beta2, amu, \
            odd = par[i:i + self.npars]

        self.alpha = par[10] - delp1 if i != 10 else par[10]
        self.coef = lambda process: self.getcoef(process, odd)
        return self.npars

    def f(self, t):
        res = self.g1 * exp(- self.beta1 * t) + \
            self.g2 * exp(- self.beta2 * t) + \
            self.g3 * (1 + t / self.t0p) ** 4
        return res ** 2

    def partial_amplitude(self, s, t, process):
        g = (-1j * s) ** (self.alpha - self.alphap * t)
        return self.coef(process) * g * self.f(t)


class SecondPomeron(FirstPomeron):
    npars = 11

    def __init__(self):
        super(PoleNumeric, self).__init__()

    def setup(self, par, i):
        self.alpha, self.alp2p, self.t0p, \
            self.gp21, self.gp22, self.gp23, \
            self.betap21, self.betap22, self.tp2, mu, \
            odd = par[i:i + self.npars]

        self.coef = lambda process: self.getcoef(process, odd)
        return self.npars

    def f(self, t):
        power = (t ** 2 + self.tp2 ** 2) ** 0.5
        res = self.gp21 * exp(- self.betap21 * t) + \
            self.gp22 * exp(- self.betap22 * power - self.tp2) + \
            self.gp23 * (1 + t / self.t0p) ** 4
        return res ** 2

    def partial_amplitude(self, s, t, process):
        g = (-1j * s) ** (self.alpha - self.alp2p * t)
        return self.coef(process) * g * self.f(t)


class Odderon(FirstPomeron):
    npars = 12

    def __init__(self):
        super(PoleNumeric, self).__init__()

    def setup(self, par, i):
        delo1, self.alp2p, self.t0p,\
            self.gp1, self.gp2, self.gp3,\
            self.betap1, self.betap2, self.tp2, self.to2, mu, \
            odd = par[i:i + self.npars]

        self.alpha = par[10] - delo1 if i != 10 else par[10]
        self.coef = lambda process: self.getcoef(process, odd)
        return self.npars

    def f(self, t):
        to = t ** 2 + self.to2 ** 2
        tp = t ** 2 + self.tp2 ** 2
        res = (
            self.gp1 * exp(- self.betap1 * to ** 0.5 - self.to2) +
            self.gp2 * exp(- self.betap2 * tp ** 0.5 - self.tp2) +
            self.gp3 * (1 + t / self.t0p) ** 4
        )
        return res ** 2

    def partial_amplitude(self, s, t, process):
        g = (-1j * s) ** (self.alpha - self.alp2p * t)
        return self.coef(process) * g * self.f(t)
