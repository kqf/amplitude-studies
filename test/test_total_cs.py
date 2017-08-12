#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter


class TestTotalCrossSection(unittest.TestCase):


    def setUp(self):
        self.cscode =  110
        self.parameters = Parameter.parameters('parameters.dat')
        data = DataPoint.read_data('alldata_v1_4.dat', [self.cscode])
        self.arguments = [(x.energy, x.t) for x in data[self.cscode]]

        with open('config/tot_cs_analytic.json') as f:
            self.conf = json.load(f)


    def checkModel(self, name, data):
        model = Eikonal(name)
        actual = [model(ss ** 2, t, self.cscode, self.parameters) for ss, t in data]
        for a, b, energy in zip(self.conf[name], actual, data):
            mymsg = 'At energy {0} GeV, the values differ, nominal: {1}, actual {2}'.format(energy, a, b)
            self.assertAlmostEqual(a, b, msg=mymsg)



    def testAnalytic(self):
        self.checkModel('analytic', self.arguments)


    def testNumeric(self):
        self.checkModel('numeric', self.arguments[0:10])


    def testCompareBoth(self):
        model1, model2 = map(Eikonal, ['analytic', 'numeric'])
        for energy, t in self.arguments[0:10]:
            a, b = model1(energy ** 2, t, self.cscode, self.parameters), model2(energy ** 2, t, self.cscode, self.parameters)
            mymsg = 'At energy {0} GeV, the values differ, analytic: {1}, numeric {2}'.format(energy, a, b)
            self.assertAlmostEqual(a, b, msg=mymsg)



