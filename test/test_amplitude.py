import unittest

from amplitude.parameter import Parameter
from amplitude.amplitude import NumericAmplitude


class TestAmplutide(unittest.TestCase):


    def test_numeric_amplitude_interface(self):
        amplitude = NumericAmplitude()
        print amplitude.npars()
