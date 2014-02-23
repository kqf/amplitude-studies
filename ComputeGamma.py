#!/usr/bin/python2.7

from ROOT import *
from DataPoint import *

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
        self.canvas.Divide(3, 2)

        [ [self.canvas.cd(i + 1),
            gPad.SetLogy(),
            gPad.SetLogx(),
            g.Draw("AP")] for i, g in enumerate(self.graphs)]
        self.canvas.Update()



def main():
    m = PerformMinimization()
    m.plot()
    raw_input('press any key ...')

if __name__ == "__main__":
    main()
else:
    print 'Module loaded'
