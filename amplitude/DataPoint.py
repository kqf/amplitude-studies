#!/usr/bin/python2

class DataPoint(object):
    def __init__(self, energy, t, observable, error, dtype):
        self.energy = energy
        self.t = t
        self.observable = observable
        self.error = error
        self.dtype = dtype

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        return 'Process %d, energy %f\n' % (self.dtype, self.energy)

