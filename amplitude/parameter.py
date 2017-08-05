#!/usr/bin/python2.7

import operator

class Parameter(object):

    def __init__(self, string):
        splitted = string.split()
        self.name = splitted[0]
        self.value, self.step, self.start, self.stop = map(float, splitted[1:])

    @staticmethod
    def read_parameters(infile):
        with open(infile) as f:
            data = map(Parameter, f)
        return data

    @staticmethod
    def parameters(infile = 'new_parameters.dat'):
        data = Parameter.read_parameters(infile)
        return map(operator.attrgetter('value'), data)
