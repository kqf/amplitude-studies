from amplitude.amplitude import NumericAmplitude


def test_numeric_amplitude_interface():
    amplitude = NumericAmplitude()
    print amplitude.npars()
