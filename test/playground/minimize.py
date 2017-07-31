#!/usr/bin/python


from scipy.optimize import minimize, rosen
import operator
from numpy.random import rand
import seaborn

import numpy as np
import matplotlib.pyplot as plt

from minuitmethod import MinuitMethod

def iterations(x, a):
	res = minimize(rosen, x, method=a, tol=1e-6)
	try:
		if res.success:
			return res.nit
		return 0
	except AttributeError:
		return 999999999

def check(data, a):
	res = np.apply_along_axis(lambda x: iterations(x, a), 0, data)
	return res.mean(), np.std(res)


def scan(n = 4, ntrials = 10):
	data = rand(n, ntrials)
	# Need jacobian:
	# , 'Newton-CG'
	# , 'dogleg'
	# , 'trust-ncg'	
	# , 'COBYLA'
	algorithms = 'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'SLSQP', MinuitMethod()

	results = [check(data, a) for a in algorithms]
	mean_iterations = {a: i for a, (i, e) in zip(algorithms, results) if i > 0.0}

	for a, niter in mean_iterations.iteritems():
		string = a, niter, n 
		print "The algorithm {0} requires {1} iterations for {2} parameters".format(*string)

	print 'The best algorithm is'
	print min(mean_iterations.iteritems(), key=operator.itemgetter(1))

	y, err = zip(*results)
	x = np.ogrid[0:len(y)]

	plt.errorbar(x, y, yerr=err)
	plt.xticks(x, mean_iterations.keys(), rotation='vertical')
	plt.yscale("log", nonposx='clip')
	plt.show()


def main():
	# for i in range(10, 40): scan(i)
	scan(40)


if __name__ == '__main__':
	main()