import pytest
import numpy as np


@pytest.mark.onlylocal
@pytest.mark.parameterize("filename, ofilename", [
    ("first.pomeron.dat", "first"),
    ("second.pomeron.dat", "second"),
])
def draw_differenence(filename, ofilename):
    import seaborn  # noqa
    import matplotlib.pyplot as plt
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
