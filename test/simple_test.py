#!/usr/bin/python

import ROOT
import json
import unittest
from amplitude.DataPoint import DataPoint
from amplitude.models import Model
from amplitude.parameter import Parameter 
import numpy as np
from scipy import integrate
from cmath import log, exp, pi
from cmath import log as cdlog
from cmath import exp as cdexp


class TestSimple(unittest.TestCase):


    def setUp(self):
        self.cscode =  110
        self.parameters = Parameter.parameters('new_parameters.dat')

        with open('config/tot_cs_analytic.json') as f:
            self.conf = json.load(f)


    @unittest.skip('')
    def testPole(self):
        s = 5.009600000000000 
        model = Model('triples_test')
        nominal = complex(934.4414471233758, -55.09380315038144)

        model.set_parameters(self.cscode, self.parameters)
        value = model.A_amplitude(s ** 2)

        print 'Nominal {0} and actual {1}'.format(nominal, value)


    @unittest.skip('')
    def testPaperFormulaForTripleExponent(self):
        # Check the formula as it was given by the author
        #
        def testAmplitude(s):
                def H(x):
                    delp1   =  0.57192E-01  
                    alp1p   =  0.38177E+00  
                    gp1     =  0.25548E+02  
                    betap11 =  0.76200E-11  
                    betap12 =  0.13120E+02  
                    betap13 =  0.24283E+01  
                    cp11    =  0.88460E+01  
                    cp12    =  0.60626E+01  
                    alp2    =  0.11609E+01  
                    alambda =  0.100000E+01


                    alp1 =  alp2 - delp1
                    cp13 = gp1

                    aim = 1j 
                    cdllp = cdlog(- aim * s)
 
                    r1pom11 = betap11 + alp1p * cdllp 
                    r1pom22 = betap12 + alp1p * cdllp 
                    r1pom33 = betap13 + alp1p * cdllp 

                    fbp1 = ( cp11 * cdexp(-0.25 * x ** 2 / r1pom11) * (0.5 / r1pom11) + 
                             cp12 * cdexp(-0.25 * x ** 2 / r1pom22) * (0.5 / r1pom22) + 
                             cp13 * cdexp(-0.25 * x ** 2 / r1pom33) * (0.5 / r1pom33) )

                    h = -fbp1*(-aim*s) ** alp1 / (8. * pi * s)


                    L = 2j * alambda
                    H = (exp(L * h) - 1) / L
                    print 'H: {0} h: {1} \t s: {2} \t b: {3} fbp1 {4}'.format(H, h, s, x, fbp1,)
                    return H

                b = 0.125000E+02
                print H(b)

                # f = lambda b: H(b).imag * b
                # imag =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

                # f = lambda b: H(b).real * b
                # real =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 
                return 0 #complex(real, imag) * (8 * pi * s)

        s = 5.009600000000000 
        nominal = complex(-55.09380315038144, 934.4414471233758)

        value = testAmplitude(s ** 2)
        print 'Nominal {0} and actual {1}'.format(nominal, value)

    def testPaperFormulaForRegularResidue(self):
        # Check the formula as it was given by the author
        #
        def testAmplitude(s):
                def H(x):
                    alf =  0.690000E+00 
                    alfp =  0.840000E+00
                    gf =  0.264894E+03 
                    bf =  0.663096E+00
                    alambda =  0.100000E+01

                    signature = -1
                    cdllp = log(-1j * s)

                    rf2 = bf + alfp * cdllp 
                    h = signature * gf * exp(alf * cdllp) * (0.5 / rf2) * exp(-0.25 * x ** 2 / rf2) / (8 * pi * s)

                    L = 2j * alambda
                    H = (exp(L * h) - 1) / L
                    print 'H: {0} h: {1} \t s: {2} \t b: {3}'.format(H, h, s, x)
                    return H

                b = 0
                # b = 0.125000E+02
                print H(b)

                # f = lambda b: H(b).imag * b
                # imag =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

                # f = lambda b: H(b).real * b
                # real =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 
                return 0
                # return complex(real, imag) * (8 * pi * s)

        s = 5.009600000000000 
        nominal = complex(-55.09380315038144, 934.4414471233758)

        value = testAmplitude(s ** 2)
        print 'Nominal {0} and actual {1}'.format(nominal, value)

