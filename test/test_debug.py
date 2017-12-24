from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter


# NB: Use this test to debug, and print intermediate results


def test_model_debug():
    model, parameters = Eikonal('triples'), Parameter.parameters(
        'triple_exp_parameters.dat')
    data = {110: [DataPoint(0.194180E+02, 0, 0, 0, 110)]}
    for d, points in data.iteritems():
        print type(points)
        for p in points:
            y = model(p.energy ** 2, p.t, d, parameters)
            print 'The result', p.energy, p.t, y
