#!/usr/bin/python

import ROOT
import unittest
from amplitude.DataPoint import DataPoint
from amplitude.amplitude import DataManager


class Draw(unittest.TestCase):

    def setUp(self):
        datasets = [110]#, 210, 310, 111, 211, 311]
        names = [ "#sigma_{pp}", "#rho_{pp}", "d#sigma_{pp}/dt", "#sigma_{ p#bar{p} }", "#rho_{ p#bar{p} }" , "d#sigma_{ p#bar{p} }/dt"]
        filename = 'alldata_v1_4.dat'
        self.names = {d: n for n, d in zip(names, datasets)}
        self.data = DataPoint.read_data(filename, datasets)
        self.approximation = DataManager(filename).get_approximation()
        
        
    def create_graph(self, datacode, data):
        """Creates TGraphErrors, with differential cross section data"""
        graph, name = ROOT.TGraphErrors(), self.names[datacode]
        graph.SetName(name)
        graph.SetTitle(name)

        for i, p in enumerate(data):
            x =  p.t if (data[0].dtype // 300 == 1) else p.energy
            y = p.observable
            if data[0].dtype // 300 == 1:
                y *= p.energy

            graph.SetPoint(i, x, y)
            graph.SetPointError(i, 0, p.error)

        graph.GetXaxis().SetTitle('#sqrt{s}, GeV')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph


    def testDraw(self):
        canvas = ROOT.TCanvas('canvas', 'Non-linear trajectories', 800, 600)

        if len(self.data) > 1:
            canvas.Divide(2, 3)

        graphs = [self.create_graph(code, data) for code, data in self.data.iteritems()]

        for i, g in enumerate(graphs):
            pad = canvas.cd(i + 1)
            g.Draw('AP')
            self.approximation.Draw('same')
            self.update_pad(pad)

        canvas.Update()
        raw_input('Press ENTER ...')

    def update_pad(self, pad):
        pad.SetTickx()
        pad.SetTicky()
        pad.SetLogy()
        pad.SetLogx()
        pad.SetGridx()
        pad.SetGridy()



