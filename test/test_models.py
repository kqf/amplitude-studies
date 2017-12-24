import pytest
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal, UMatrix
from amplitude.parameter import Parameter


MODELS = {
    'analytic': 'parameters.dat',
    'numeric': 'parameters.dat',
    'triples': 'new_parameters.dat'
}


@pytest.fixture(scope="module")
def models():
    return {
        modelname: Parameter.parameters(filename)
        for modelname, filename in MODELS.iteritems()
    }


@pytest.fixture(scope="module")
def data():
    return DataPoint.read_data("dsdtout.dat")


@pytest.mark.parametrize("unitarization, utype", [
    ('eikonal', Eikonal),
    ('u-matrix', UMatrix),
])
def dump(unitarization, utype, models, data, npoints=2):
    for m, par in models.iteritems():
        print
        print 'Calculating {0}:'.format(m)
        model = utype(m)
        for d, points in data.iteritems():
            for p in points[0:npoints]:
                y = model(p.energy ** 2, p.t, d, par)
                msg = 'For {} {} model observable = {}, theoretical value = {}'
                print msg.format(m, unitarization, p.observable, y)
