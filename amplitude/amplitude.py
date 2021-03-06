
from poles import Pole, PoleNumeric, TripleExponentPole
import poles


class NumericAmplitude(object):
    def __init__(self):
        super(NumericAmplitude, self).__init__()
        self.poles = PoleNumeric(), PoleNumeric(), PoleNumeric()
        self.type = PoleNumeric

    def h(self, s, x, code):
        return self.type.h_amplitude(self.poles, s, x, code)

    def set_parameters(self, par):
        i, par = 0, [p for p in par]  # ROOT provides buffer, not array
        for poleNumeric in self.poles:
            i += poleNumeric.setup(par, i)
        llambda = par[-1]
        return llambda

    def use_single_pole(self, i):
        self.poles = self.poles[i:i + 1]

    def npars(self):
        return sum(i.npars for i in self.poles)


class AnalyticAmplitude(NumericAmplitude):
    def __init__(self):
        super(AnalyticAmplitude, self).__init__()
        self.poles = Pole(), Pole(), Pole()
        self.type = Pole


class TripleExponentAmplitude(NumericAmplitude):
    def __init__(self):
        super(TripleExponentAmplitude, self).__init__()
        self.poles = TripleExponentPole(), TripleExponentPole(
        ), TripleExponentPole(), Pole(), Pole()
        self.type = Pole


class NonlinearAmplitude(NumericAmplitude):
    def __init__(self):
        super(NonlinearAmplitude, self).__init__()
        self.poles = poles.FirstPomeron(), poles.SecondPomeron(
        ), poles.Odderon(), PoleNumeric(), PoleNumeric()
        self.type = PoleNumeric
