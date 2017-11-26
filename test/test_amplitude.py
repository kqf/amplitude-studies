import unittest

from amplitude.amplitude import NumericAmplitude
from amplitude.parameter import Parameter 


class TestAmplutide(unittest.TestCase):


	def test_numeric_amplitude_interface(self):
		amplitude = NumericAmplitude()
		print amplitude.npars()

		self.amplitude()




