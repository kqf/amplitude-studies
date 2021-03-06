import operator


class Parameter(object):

    def __init__(self, string):
        splitted = string.split()
        self.name = splitted[0]
        self.value, self.step, self.start, self.stop = map(float, splitted[1:])

    @staticmethod
    def read_parameters(infile):
        # Store all the data in one folder, fore this behaviour
        with open('input-data/' + infile) as f:
            data = map(Parameter, f)
        return data

    @staticmethod
    def parameters(infile='new_parameters.dat'):
        data = Parameter.read_parameters(infile)
        return map(operator.attrgetter('value'), data)
