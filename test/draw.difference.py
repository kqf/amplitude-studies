#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import seaborn


# TODO: rewrite as test
#

def draw_differenence(filename, ofilename):
    e, b, re, im = np.loadtxt('test/' + filename).T
    bb, rre, iim = np.loadtxt('test/m.' + filename).T

    plt.plot(b, (re - rre), 'o', label='real')
    plt.plot(b, (im - iim), 'o', label='imag')
    plt.title(filename.replace('.', ' '))
    plt.ylabel('h1 - h2')
    plt.xlabel('impact parameter, b')
    plt.legend()
    plt.savefig(ofilename + '.png')
    plt.show()


def main():
    draw_differenence("first.pomeron.dat", "first")
    draw_differenence("second.pomeron.dat", "second")


if __name__ == '__main__':
    main()
