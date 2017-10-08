#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter



# NB: Use this test to debug, and print intermediate results


class TestIntermediateResults(unittest.TestCase):


    def test_linear_model(self):
        test_variables = 110,
        model, parameters = Eikonal('triples'), Parameter.parameters('triple_exp_parameters.dat')
        data = DataPoint.read_data('dsdtout.dat')
        # inp_data = {tv: data[tv] for tv in data if tv in test_variables}
        inp_data = {110: [DataPoint(0.194180E+02, 0, 0, 0, 110)] }
        self.run_model(model, parameters, inp_data)


    def run_model(self, model, parameters, data):
        for d, points in data.iteritems():
            print type(points)
            for p in points:
                y = model(p.energy ** 2, p.t, d, parameters)
                print 'The result', p.energy, p.t, y

