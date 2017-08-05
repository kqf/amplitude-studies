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


    def testPole(self):
        s = 5.009600000000000 
        model = Model('triples_test')
        nominal = complex(934.4414471233758, -55.09380315038144)

        model.set_parameters(self.cscode, self.parameters)
        value = model.A_amplitude(s ** 2)

        print 'Nominal {0} and actual {1}'.format(nominal, value)


    def testPaperFormula(self):
        # Check the formula as it was given by the author
        #

        def testAmplitude(s):
                def H(x):
                    delp1   = 0.571918E-01        
                    alp1p   = 0.381769E+00        
                    gp1     = 0.255478E+02        
                    betap11 = 0.762002E-11        
                    betap12 = 0.131201E+02        
                    betap13 = 0.242827E+01        
                    cp11    = 0.884597E+01        
                    cp12    = 0.606263E+01        
                    alp2    = 0.116093E+01        
                    alambda = 0.100000E+01

                    alp1 = alp2 - delp1
                    cp13 = gp1

                    aim = 1j 
                    cdllp = cdlog(-aim * s)
 
                    r1pom11 = betap11 + alp1p * cdllp 
                    r1pom22 = betap12 + alp1p * cdllp 
                    r1pom33 = betap13 + alp1p * cdllp 

                    fbp1 = ( cp11 * cdexp(-0.25 * x ** 2 / r1pom11) * (0.5 / r1pom11) + 
                             cp12 * cdexp(-0.25 * x ** 2 / r1pom22) * (0.5 / r1pom22) + 
                             cp13 * cdexp(-0.25 * x ** 2 / r1pom33) * (0.5 / r1pom33) )

                    h = -fbp1*(-aim*s) ** alp1 / (8 * pi * s)



                    # ss = -1j * s
                    # alphalog = log(ss) * alp1p

                    # f = lambda r: r + alphalog 
                    # r1, r2, r3 = map(f, (betap11, betap12, betap13))

                    # fbp1 =(cp11 * exp(-0.25 * x **2 / r1) * (0.5 / r1) + 
                    #        cp12 * exp(-0.25 * x **2 / r2) * (0.5 / r2) + 
                    #        cp13 * exp(-0.25 * x **2 / r3) * (0.5 / r3))

                    # h = -fbp1 * ss ** alp1 / (8 * pi * s)
                    L = 2j * alambda
                    H = (exp(L * h) - 1) / L
                    return H

                f = lambda b: H(b).imag * b
                imag =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 

                f = lambda b: H(b).real * b
                real =  integrate.quad(f, 0, 40)[0]  # integral  zero to lower 
                return complex(real, imag) * (8 * pi * s)

        s = 5.009600000000000 
        nominal = complex(934.4414471233758, -55.09380315038144)

        value = testAmplitude(s ** 2)
        print 'Nominal {0} and actual {1}'.format(nominal, value)

