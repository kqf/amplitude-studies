
from poles import Pole, PoleNumeric, TripleExponentPole

class NumericAmplitude(object):
    def __init__(self):
        super(NumericAmplitude, self).__init__()
        self.poles = PoleNumeric(), PoleNumeric(), PoleNumeric()
        self.type = PoleNumeric

    def h(self, s, x):
        return self.type.h_amplitude(self.poles, s, x)

    def set_parameters(self, code, par):
        i, par = 0, [p for p in par] # ROOT provides buffer, not array
        for poleNumeric in self.poles:
            i += poleNumeric.setup(par, i, code)
        llambda = par[-1]
        return llambda


class AnalyticAmplitude(NumericAmplitude):
    def __init__(self):
        super(AnalyticAmplitude, self).__init__()
        self.poles = Pole(), Pole(), Pole()
        self.type = Pole
                        

class TripleExponentAmplitude(NumericAmplitude):
    def __init__(self):
        super(TripleExponentAmplitude, self).__init__()
        self.poles = TripleExponentPole(), TripleExponentPole(), TripleExponentPole(), Pole(), Pole()
        self.type = Pole

