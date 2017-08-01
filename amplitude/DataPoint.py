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
    def read_dataset(infile, dataset):
        """Reading data from alldata_v1_4.dat file
           0 -- \sqrt(s)
           1 -- -t         
           2 -- observable 
           5 -- error
           6 -- id of the observable
        """
        raw_data = []
        with open(infile,'r') as file:
            for line in file:
                data = line.lower().split()

                dataset_in_file = int( float(data[6]))

                if dataset_in_file != dataset:
                    continue

                energy, t, observable, error = map(float, data[0:4])

                raw_data.append( DataPoint(energy, t,
                    observable, error, dataset) )

        return raw_data

    @staticmethod
    def read_data(infile, datasets):
        raw_data = {d: DataPoint.read_dataset(infile, d) for d in datasets}
        return raw_data



