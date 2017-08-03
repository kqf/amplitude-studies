#!/usr/bin/python2.7

from ROOT import *
import ROOT
from DataPoint import DataPoint
from models import Model


class DataManager(object):
    datasets = [110, 210, 310, 111, 211, 311]

    def __init__(self, infile):
        self.data = DataPoint.read_data(infile, self.datasets)
        self.model = Model()


    def get_approximation(self):
        # TODO: write normal input data handler
        parameters = [40.3043, 1.10517, 0.35, 2.32537, 1, 
                      117.221, 0.791348, 1.31638, 1.98679, 1,
                      102.76, 0.5, 1.2, 8.80651, 1, 
                      0.5]

        # self.model([4], parameters)
        # TODO: remove hardocded numbers
        func = ROOT.TF1('func', lambda x, p: self.model(x[0], 110, p), 5, 3000, len(parameters))

        for i, p in enumerate(parameters): 
            func.SetParameter(i, p)

        return func

