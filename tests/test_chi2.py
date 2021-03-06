import pytest
from amplitude.datapoint import DataPoint
from amplitude.models import Eikonal
from amplitude.parameter import Parameter


def accept_point(p, datacode):

    if datacode // 300 == 0:
        return 5.000 < p.energy and p.energy < 8000.450

    return (
        (19.0 < p.energy and p.energy < 8000.0) and
        (0.010 < p.t and p.t < 15.000)
    )


def run_chi2(model, parameters, data):
    chi2, npoints = 0, 0
    for d, points in data.iteritems():
        accepted = (p for p in points if accept_point(p, d))
        chi2_, npoints_ = 0, 0
        for p in accepted:
            y = model(p.energy ** 2, p.t, d, parameters)

            chi2__ = ((p.observable - y) / p.error) ** 2
            # print p.energy, p.observable, y, p.error,  chi2__

            chi2_ += chi2__
            # print p.energy, p.t,  chi2__, d
            npoints_ += 1
        chi2 += chi2_
        npoints += npoints_
        msg = 'Chi^2\t{0} chi^2/ndf \t{3} for " \
            "\t{1} \tpoints, per process {2}'
        print msg.format(chi2_, npoints_, d, chi2_ / npoints_)

    print
    print 'Chi^2     {0} for {1} points'.format(chi2, npoints)
    print 'Chi^2/ndf {0} for {1} points'.format(chi2 / npoints, npoints)


@pytest.mark.onlylocal
def test_linear_model():
    model, parameters = Eikonal('triples'), Parameter.parameters(
        'triple_exp_parameters.dat')
    data = DataPoint.read_data('dsdtout.dat')
    run_chi2(model, parameters, data)
