#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.DataPoint import DataPoint
from amplitude.models import Model


class TestTotalCrossSection(unittest.TestCase):


    def setUp(self):
        self.cscode =  110
        self.parameters = [40.3043, 1.10517, 0.35, 2.32537, 1, 117.221, 0.791348, 1.31638, 1.98679, 1, 102.76, 0.5, 1.2, 8.80651, -1, 0.5]
        data = DataPoint.read_data('alldata_v1_4.dat', [self.cscode])
        self.energy = [x.energy for x in data[self.cscode]]

        with open('config/tot_cs_analytic.json') as f:
            self.conf = json.load(f)


    def checkModel(self, name, data):
        model = Model(name)
        actual = [model(x, self.cscode, self.parameters) for x in data]
        for a, b, energy in zip(self.conf[name], actual, data):
            mymsg = 'At energy {0} GeV, the values differ, nominal: {1}, actual {2}'.format(energy, a, b)
            self.assertAlmostEqual(a, b, msg=mymsg)



    def testAnalytic(self):
        self.checkModel('analytic', self.energy)


    def testNumeric(self):
        self.checkModel('numeric', self.energy[0:10])


    @unittest.skip('')
    def testCompareBoth(self):
        model1, model2 = map(Model, ['analytic', 'numeric'])
        for energy in self.energy:
            a, b = model1(energy, self.cscode, self.parameters), model2([energy], self.cscode, self.parameters)
            mymsg = 'At energy {0} GeV, the values differ, analytic: {1}, numeric {2}'.format(energy, a, b)
            self.assertAlmostEqual(a, b, msg=mymsg)



