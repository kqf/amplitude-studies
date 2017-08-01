#!/usr/bin/python2.7

from ROOT import *
import ROOT
from DataPoint import DataPoint
from Formulas import fit_function


class DataManager(object):
    datasets = [110, 210, 310, 111, 211, 311]

    def __init__(self, infile):
        self.data = DataPoint.read_data(infile, self.datasets)


    def plot_approximation(self):
        canvas = ROOT.TCanvas('canvas', 'Non-linear trajectories', 800, 600)
        # TODO: write normal input data handler
        parameters = [40.3043, 117.221, 102.76,0.35,1.10517,1.31638,0.791348,1.2,0.5,0.160351,2.32537,1.98679,8.80651]

        # TODO: remove hardocded numbers
        func = ROOT.TF1('func', fit_function, 5, 3000, len(parameters))
        [ func.SetParameter(i, p) for i, p in enumerate(parameters) ]

        func.Draw('APL')
        canvas.cd(1)
        canvas.Update()



def main():
    m = PerformMinimization()
    m.plot_approximation()

if __name__ == "__main__":
    main()
