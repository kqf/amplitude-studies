#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.DataPoint import DataPoint
from amplitude.models import Eikonal, UMatrix
from amplitude.parameter import Parameter


class TestTotalCrossSection(unittest.TestCase):


    def setUp(self):
        models = {'analytic': 'parameters.dat', 'numeric': 'parameters.dat', 'triples': 'new_parameters.dat'}

        self.models = {m: Parameter.parameters(p) for m, p in models.iteritems()}
        self.data = DataPoint.read_data('dsdtout.dat')
        self.npoints = 2

    def dump(self, unitarization, utype):
        for m, par in self.models.iteritems():
            print
            print 'Calculating {0}:'.format(m)
            model = utype(m)
            for d, points in self.data.iteritems():
                for p in points[0:self.npoints]:
                    y = model(p.energy ** 2, p.t, d, par) 
                    print 'For {0} {1} model observable = {2}, theoretical value = {3}'.format(m, 
                        unitarization, p.observable, y)


    def testEikonal(self):
        self.dump('eikonal', Eikonal)


    def testUMatrix(self):
        self.dump('u-matrix',UMatrix)
