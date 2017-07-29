#!/usr/bin/python

import ROOT
import numpy as np
from scipy.optimize import minimize, rosen, OptimizeResult
from array import array


class MinuitMethod(object):

	def __init__(self, xtol = 1e-6):
		super(MinuitMethod, self).__init__()
		self.xtol = xtol
		self.nit = 0


	def fcn(self, npar, gin, f, par, iflag):
		# This is a nuissance of ROOT interface
		self.nit += 1
		pars = [par[i] for i in range(self.npars)]
		return self.func(pars)


	def __call__(self, fun, x0, args, **kwargs):
		self.func = fun
		self.npars = len(x0)
		minuit = ROOT.TMinuit(self.npars)
		minuit.SetFCN(self.fcn)
		ierflg = ROOT.Long(2999)
		for i, p in enumerate(x0):
			minuit.mnparm(i, str(p), p, self.xtol, 0, 0, ierflg)

		# Setup the minuit parameters
		arglist = array('d', 10*[0.])

		arglist[0] = 1
		minuit.mnexcm("SET ERR", arglist, 1, ierflg)

		arglist[0] = 500
		arglist[1] = 1.
		self.nit = 0
		minuit.mnexcm("MIGRAD", arglist, 2, ierflg)
		# result = OptimizeResult(fun=fval, direc=direc, nit=iter, nfev=fcalls[0], status=warnflag, success=(warnflag == 0), message=msg, x=x)

		amin, edm, errdef = ROOT.Double(0.18), ROOT.Double(0.19), ROOT.Double(0.20)
		nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
		minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat)

		success = amin < 2.0 and edm < kwargs['tol']
		msg = 'Success' if success else 'Fail'

		optpar = [0] * self.npars
		opterr = [0] * self.npars

		# This will not work due to pyROOT limitations
		# for i, (x, ex) in enumerate(zip(optpar, opterr)):

		for i in range(self.npars):
			x, ex = ROOT.Double(0), ROOT.Double(0)
			minuit.GetParameter(i, x, ex)
			optpar[i] = x
			opterr[i] = ex

		return OptimizeResult(nit=self.nit, success=success, message = msg, x=optpar)

	def __str__(self):
		return 'Minuit'


def main():
	x = np.random.rand(4)
	res = minimize(rosen, x, method=MinuitMethod(), tol=1e-6)
	
	print MinuitMethod()

if __name__ == '__main__':
	main()