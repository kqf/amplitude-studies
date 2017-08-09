#!/usr/bin/python

import ROOT
import unittest
from amplitude.datapoint import DataPoint
from amplitude.parameter import Parameter
from amplitude.models import Eikonal


class Draw(unittest.TestCase):

    def setUp(self):
        # datasets = [210] # [110, 111, 210, 211]#, 310, 311]
        datasets =  [110, 111, 210, 211, 310, 311][::-1]
        names = [ "#sigma_{pp}", "#sigma_{ p#bar{p} }", "#rho_{pp}", "#rho_{ p#bar{p} }", "d#sigma_{pp}/dt" , "d#sigma_{ p#bar{p} }/dt"]
        self.filename = 'alldata_v1_4.dat'
        self.names = {d: n for n, d in zip(names, datasets)}
        self.data = DataPoint.read_data(self.filename, datasets)
        self.model = Eikonal('triples')


    def approximation(self, code, energy = 19.4):
        par = Parameter.parameters('new_parameters.dat')
        if code // 300 == 0:
            f = lambda x, p: self.model(x[0] ** 2, 0, code, p)
            func = ROOT.TF1('func', f, 5, 1e5, len(par)) 
        else:
            f = lambda x, p: self.model(energy ** 2, x[0], code, p)
            func = ROOT.TF1('func', f, 0, 15, len(par))

        for i, p in enumerate(par): 
            func.SetParameter(i, p)
        return func 
        
    def graph_vs_approx(self, datacode, data):
        """Creates TGraphErrors, with differential cross section data"""
        graph, name = ROOT.TGraphErrors(), self.names[datacode]
        graph.SetName(name)
        graph.SetTitle(name)

        for i, p in enumerate(data):
            x =  p.t if (data[0].dtype // 300 == 1) else p.energy
            y = p.observable
            if data[0].dtype // 300 == 1:
                if p.energy > 19.2 and p.energy < 19.4 :
                    continue
                # y *= p.energy

            graph.SetPoint(i, x, y)
            graph.SetPointError(i, 0, p.error)

        graph.GetXaxis().SetTitle('#sqrt{s}, GeV')
        graph.GetYaxis().SetTitle('#frac{d#sigma}{dt}, mb/GeV^{2}')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(46)
        return graph, self.approximation(datacode)


    def testDraw(self):
        canvas = ROOT.TCanvas('canvas', 'Non-linear trajectories', 800, 600)

        for code in sorted(self.data, reverse = True):
            data = self.data[code]
            graph, approx = self.graph_vs_approx(code, data)
            pad = canvas.cd()
            self.update_pad(pad, code)
            graph.Draw('AP')
            approx.Draw('same')
            canvas.Update()
            canvas.SaveAs(str(code) + '.pdf')
            raw_input('Press ENTER ...')


    def update_pad(self, pad, code):
        pad.SetLogy(code // 300 == 1)
        pad.SetLogx(code // 300 == 0)

        pad.SetTickx()
        pad.SetTicky()
        # pad.SetLogy()
        pad.SetLogx()
        pad.SetGridx()
        pad.SetGridy()
