import pytest
import json
import numpy as np
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter


@pytest.fixture(scope="module")
def nominal():
    with open('config/tot_cs_analytic.json') as f:
        conf = json.load(f)
    return conf


def data(cscode=110, datatype="analytic"):
    data = DataPoint.read_data('alldata_v1_4.dat', [cscode])
    dataset = [(x.energy, x.t) for x in data[cscode]]
    if datatype == "numeric":
        return dataset[0:10]
    return dataset


@pytest.fixture(scope="module")
def parameters(infile="parameters.dat"):
    return Parameter.parameters('parameters.dat')


@pytest.mark.skip("")
@pytest.mark.parametrize("name", [
    'analytic'
    'numeric',
])
def test_calculate_model(name, parameters, nominal, cscode=110):
    model, indata = Eikonal(name), data(name)
    actual = [model(ss ** 2, t, cscode, parameters)
              for ss, t in indata]
    for a, b, energy in zip(nominal[name], actual, indata):
        msg = 'At energy {} GeV, the difference, nominal: {}, actual {}'
        assert pytest.approx(a) == b, msg.format(energy, a, b)


@pytest.mark.onlylocal
def test_models_are_consistent(parameters, cscode=110):
    model1, model2 = map(Eikonal, ['analytic', 'numeric'])
    for energy, t in data("numeric"):
        a = model1(energy ** 2, t, cscode, parameters),
        b = model2(energy ** 2, t, cscode, parameters)
        msg = 'At energy {} GeV, difference analytic: {}, numeric {}'
        np.testing.assert_almost_equal(a, b, err_msg=msg.format(energy, a, b))
