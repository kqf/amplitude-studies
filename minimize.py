#!/usr/bin/python


from scipy.optimize import minimize, rosen
import operator
from numpy.random import rand


def scan(n = 4):
	x0 = rand(n)
	# x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
	# Need jacobian:
	# , 'Newton-CG'
	# , 'dogleg'
	# , 'trust-ncg'	
	# , 'COBYLA'
	algorithms = 'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'SLSQP'


	f = lambda x: minimize(rosen, x0, method=x, tol=1e-6)
	data = {a: f(a) for a in algorithms}
	iterations = {a: r.nit for a, r in data.iteritems()}

	for a, niter in iterations.iteritems():
		string = a, niter, len(x0)
		print "The algorithm {0} requires {1} iterations for {2} parameters".format(*string)

	print 'The best algorithm is'
	print min(iterations.iteritems(), key=operator.itemgetter(1))


def main():
	for i in range(10, 40):
		scan(i)


if __name__ == '__main__':
	main()