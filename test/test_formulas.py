import unittest
from amplitude.models import Eikonal
from amplitude.parameter import Parameter

import numpy as np
from scipy import integrate
from cmath import log, exp, pi
from cmath import log as cdlog
from cmath import exp as cdexp


# Test the formula
#
def checkTriplePoleFormula(s):
    def H(x):
        delp1 = 0.57192E-01
        alp1p = 0.38177E+00
        gp1 = 0.25548E+02
        betap11 = 0.76200E-11
        betap12 = 0.13120E+02
        betap13 = 0.24283E+01
        cp11 = 0.88460E+01
        cp12 = 0.60626E+01
        alp2 = 0.11609E+01
        alambda = 0.100000E+01

        alp1 = alp2 - delp1
        cp13 = gp1

        aim = 1j
        cdllp = cdlog(- aim * s)

        r1pom11 = betap11 + alp1p * cdllp
        r1pom22 = betap12 + alp1p * cdllp
        r1pom33 = betap13 + alp1p * cdllp

        fbp1 = (cp11 * cdexp(-0.25 * x ** 2 / r1pom11) * (0.5 / r1pom11) +
                cp12 * cdexp(-0.25 * x ** 2 / r1pom22) * (0.5 / r1pom22) +
                cp13 * cdexp(-0.25 * x ** 2 / r1pom33) * (0.5 / r1pom33))

        h = -fbp1 * (-aim * s) ** alp1 / (8. * pi * s)

        L = 2j * alambda
        H = (exp(L * h) - 1) / L
        return H

    # b = 0.125000E+02
    # print H(b)

    def f(b):
        return H(b).imag * b
    imag = integrate.quad(f, 0, 40)[0]  # integral  zero to lower

    def f(b):
        return H(b).real * b
    real = integrate.quad(f, 0, 40)[0]  # integral  zero to lower
    # return 0
    return complex(real, imag) * (8 * pi * s)


def dump_to_file(h, filename='m.second.pomeron.dat'):
    b = np.arange(0, 5., 0.1)
    data = [h(ib)[1] for ib in b]
    real, imag = np.array([i.real for i in data]), np.array(
        [i.imag for i in data])
    np.savetxt('test/m.' + filename + '.dat', np.asarray([b, real, imag]).T)


def checkTriplePoleFormulaNew(s):
    def h(x, verbose=False):
        # Genral constants
        alp2 = 0.11581E+01
        alambda = 0.100000E+01

        # # First Pomeron
        # delp1  = 0.53465E-01
        # alp1p  = 0.37310E+00

        # gp11  = 0.91325E+01
        # gp12  = 0.57979E+01
        # gp13  = 0.25445E+02

        # betap11  = 0.10001E+00
        # betap12  = 0.13518E+02
        # betap13  = 0.25251E+01
        # alp1 =  alp2 - delp1

        # Second Pomeron
        # alp1 = 0.11581E+01

        # alp1p = 0.82485E-01
        # gp11 = 0.40327E+01
        # gp12 = 0.19212E+00
        # gp13 = 0.48001E+01

        # betap11 = 0.10499E+01
        # betap12 = 0.10000E+00
        # betap13 = 0.40484E+01

        # Odderon
        delp1 = 0.42180E-01
        alp1p = 0.55991E-01
        gp11 = 0.10738E+00
        gp12 = 0.13869E+01
        gp13 = 0.47682E-05
        betap11 = 0.10001E+00
        betap12 = 0.41992E+01
        betap13 = 0.89147E+01
        alp1 = alp2 - delp1

        aim = 1j
        cdllp = cdlog(- aim * s)

        r1pom11 = betap11 + alp1p * cdllp
        r1pom22 = betap12 + alp1p * cdllp
        r1pom33 = betap13 + alp1p * cdllp

        fbp1 = (gp11 * cdexp(-0.25 * x ** 2 / r1pom11) * (0.5 / r1pom11) +
                gp12 * cdexp(-0.25 * x ** 2 / r1pom22) * (0.5 / r1pom22) +
                gp13 * cdexp(-0.25 * x ** 2 / r1pom33) * (0.5 / r1pom33))

        h = -fbp1 * (-aim * s) ** alp1 / (8. * pi * s)

        L = 2j * alambda
        H = (exp(L * h) - 1) / L
        return H, h

    dump_to_file(h, 'odderon')

    def H(b):
        return h(b)[0]

    def f(b):
        return H(b).imag * b
    imag = integrate.quad(f, 0, 40)[0]  # integral  zero to lower

    def f(b):
        return H(b).real * b
    real = integrate.quad(f, 0, 40)[0]  # integral  zero to lower
    # return 0
    return complex(real, imag) * (8 * pi * s)


def checkRegularPoleFormula(s):
    def h(x):
        # f reggeon
        alf = 0.690000E+00
        alfp = 0.840000E+00
        gf = 0.264894E+03
        bf = 0.663096E+00
        alambda = 0.100000E+01

        signature = -1
        cdllp = log(-1j * s)

        rf2 = bf + alfp * cdllp
        h = signature * gf * exp(alf * cdllp) * (0.5 / rf2) * \
            exp(-0.25 * x ** 2 / rf2) / (8 * pi * s)

        L = 2j * alambda
        H = (exp(L * h) - 1) / L
        # print 'H: {0} h: {1} \t s: {2} \t b: {3}'.format(H, h, s, x)
        return H, h

    dump_to_file(h, 'omega')

    def H(b):
        return h(b)[0]

    def f(b):
        return H(b).imag * b
    imag = integrate.quad(f, 0, 40)[0]  # integral  zero to lower

    def f(b):
        return H(b).real * b
    real = integrate.quad(f, 0, 40)[0]  # integral  zero to lower
    # return 0
    return complex(real, imag) * (8 * pi * s)


#  NB: Don't use assertions here, as this code and values will change
#
@unittest.skip('')
def test_calculates_single_pole(cscode=110):
    s, t = 5.009600000000000, 0
    model = Eikonal('triples')
    nominal = complex(-630.265349098, 2542.98388891)

    model.set_parameters(Parameter.parameters('new_parameters.dat'))
    value = model.A_amplitude(s ** 2, t, cscode, False)
    print 'Total amplitude: Nominal {0} and actual {1}'.format(nominal, value)


@unittest.skip('')
def test_formula_from_triple_exp_paper():
    # Check the formula as it was given by the author
    #
    s = 5.009600000000000
    nominal = complex(206.494887954, 1207.64861523)

    value = checkTriplePoleFormula(s ** 2)
    print 'Triple pole: Nominal {0} and actual {1}'.format(nominal, value)


@unittest.skip('')
def test_regular_residuer_form_paper():
    # Check the formula as it was given by the author
    #

    s = 5.009600000000000
    nominal = complex(-625.921351373, 1952.41665889)

    value = checkRegularPoleFormula(s ** 2)
    print 'Regular pole:  Nominal {0} and actual {1}'.format(
        nominal, value)


def test_new_parameters():
    sqrts = 0.194180E+02
    model, parameters = Eikonal('triples'), Parameter.parameters(
        'triple_exp_parameters.dat')
    # NB: First set parameters and then delete poles
    model.set_parameters(parameters)
    model.amplitude.use_single_pole(0)
    formula = checkTriplePoleFormulaNew(sqrts ** 2)

    print "Explicit formula: {} {}".format(
        formula,
        model.A_amplitude(sqrts ** 2, 0, 110, skipreal=False))

    checkRegularPoleFormula(sqrts ** 2)
