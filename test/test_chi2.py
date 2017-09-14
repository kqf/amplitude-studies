#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter


class TestChi2(unittest.TestCase):


    def test_linear_model(self):
        model, parameters = Eikonal('triples'), Parameter.parameters('triple_exp_parameters.dat')
        data = DataPoint.read_data('dsdtout.dat')
        self.run_chi2(model, parameters, data)

    @unittest.skip('check the parameters')
    def test_nonlinear_model(self):
        model, parameters = Eikonal('nmod'), Parameter.parameters('nonlin_parameters.dat')
        data = DataPoint.read_data('dsdtout.dat')
        self.run_chi2(model, parameters, data)

    def accept_point(self, p, datacode):
        
        if datacode // 300 == 0:
            return 5.000 < p.energy and p.energy < 8000.450

        return 19.0 < p.energy and p.energy < 8000.0 and 0.010 < p.t and p.t < 15.000


    def run_chi2(self, model, parameters, data):
        chi2, npoints = 0, 0
        for d, points in data.iteritems():
            accepted = (p for p in points if self.accept_point(p, d))
            chi2_, npoints_ = 0, 0
            for p in accepted:
                y = model(p.energy ** 2, p.t, d, parameters)

                chi2__ = ((p.observable - y)/ p.error) ** 2
                # print p.energy, p.observable, y, p.error,  chi2__

                chi2_ +=  chi2__
                # print p.energy, p.t,  chi2__, d
                npoints_ += 1
            chi2 += chi2_
            npoints += npoints_

            print 'Chi^2\t{0} chi^2/ndf \t{3} for \t{1} \tpoints, per process {2}'.format(chi2_, npoints_, d, chi2_ / npoints_)


        print
        print 'Chi^2     {0} for {1} points'.format(chi2, npoints)
        print 'Chi^2/ndf {0} for {1} points'.format(chi2 / npoints, npoints)
