#!/usr/bin/python2

class DataPoint(object):
    def __init__(self, energy, t, observable, error, dtype):
        self.energy = energy
        self.t = t
        self.observable = observable
        self.error = error
        self.dtype = dtype

    def __lt__(self, other):
        return self.t < other.t

    def __repr__(self):
        return 'Process %d, energy %f\n' % (self.dtype, self.energy)

    @staticmethod
    def _extract_fields(string, filename, required):
        dataset = int(string[6])

        if dataset != required:
            return None

        if not 'dsdtout.dat' in filename:
            return map(float, string[0:4]) + [dataset]

        return map(float, string[0:2] + string[3:5]) + [dataset]

    @classmethod
    def read_dataset(klass, infile, dataset):
        f = lambda x: klass._extract_fields(x, infile, dataset)

        with open('input-data/' + infile, 'r') as ifile:
            fields = [f(line.split()) for line in ifile]

        datapoints = [klass(*i) for i in fields if i]
        return datapoints

    @classmethod
    def read_data(klass, infile, datasets = [110, 111, 210, 211, 310, 311]):
        data = {d: klass.read_dataset(infile, d) for d in datasets}
        return data



