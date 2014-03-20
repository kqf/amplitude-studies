#!/usr/bin/python2.7

from ROOT import *
from DataPoint import *
from Formulas import fit_function

class PerformMinimization(object):
    canvas = TCanvas('canvas', 'Non linear trajectories', 800, 600)
    datasets = [110, 210, 310, 111, 211, 311]

    def __init__(self):
        self.data = []
        self.title = "#sigma_{pp}"
        [ self.read_data(i) for i in self.datasets ]

        self.names = [ "#sigma_{pp}", "#rho_{pp}", "d#sigma_{pp}/dt",
                  "#sigma_{ p#bar{p} }", "#rho_{ p#bar{p} }" , "d#sigma_{ p#bar{p} }/dt"]

        self.graphs = [ self.create_graph(self.names[i], d) for i, d in enumerate(self.data)]
        self.data_is_drawn = False

    def read_data(self, dataset):
        """Reading data from alldata_v1_4.dat file
           0 -- \sqrt(s)
           1 -- -t         
           2 -- observable 
           5 -- error
           6 -- id of the observable
        """
        raw_data = []
        with open('alldata_v1_4.dat','r') as file:
            for line in file:
                data = line.lower().split()

                dataset_in_file = int( float(data[6]))

                if dataset_in_file != dataset:
                    continue

                energy = float(data[0])
                t      = float(data[1])

                if self.draw_cut(energy, t):
                    continue

                observable = float(data[2])
                error      = float(data[3])

                raw_data.append( DataPoint(energy, t,
                    observable, error, dataset) )

        self.data.append(raw_data)


    def draw_cut(self, energy, t):
        """Check if data is within draw range"""
        return False

    def create_graph(self, name, data):
        """Creates TGraphErrors, with differential cross section data"""
        # TODO: Check if data red properly
        graph = TGraphErrors()
        graph.SetName(name)
        graph.SetTitle(name + " graph")

        if data[0].dtype // 300 == 1:
            [ graph.SetPoint(i, p.t, p.observable) for i, p in enumerate(data) ]
        else:
            [ graph.SetPoint(i, p.energy, p.observable) for i, p in enumerate(data) ]

        [ graph.SetPointError(i, 0, p.error) for i, p in enumerate(data) ]

        graph.GetXaxis().SetTitle('#sqrt{s}, GeV')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)

        return graph

    def plot(self):
        # self.canvas.Divide(3, 2)
        # self.canvas.Divide(3, 2)

        self.canvas.cd(1)
        gPad.SetLogy()
        gPad.SetLogx()
        self.graphs[0].Draw("AP")



        # [ [self.canvas.cd(i + 1),
            # gPad.SetLogy(),
            # gPad.SetLogx(),
            # g.Draw("AP")] for i, g in enumerate(self.graphs)]


        self.canvas.Update()
        self.data_is_drawn = True # needed to keep track of draw parameters

    def plot_approximation(self):
        func_draw_parameter = 'same' if self.data_is_drawn else 'APL'
        self.canvas.cd(1)
        # TODO: write normal input data handler
        
        parameters = [ 40.3043,117.221,102.76,0.35,1.10517,1.31638,0.791348,1.2,0.5,0.160351,2.32537,1.98679,8.80651]

        # TODO: remove hardocded numbers
        func = TF1('func',fitFunc, 5, 3000, len(parameters))
        [ func.SetParameter(i, p) for i, p in enumerate(parameters) ]

        func.Draw('APL same' if self.data_is_drawn else 'APL')
        self.canvas.cd(1)
        self.canvas.Update()

def fitFunc(x, p):
    return fit_function(x, p)


def main():
    m = PerformMinimization()
    m.plot()
    m.plot_approximation()
    raw_input('Press any key ...')

if __name__ == "__main__":
    main()
else:
    print 'Module loaded'
