import ROOT
from ROOT import *
from array import array;

def func(x, y, par):
   value = ( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y)
   return value


class MinuitCheck(object):
   def __init__(self, *args):
      super(MinuitCheck, self).__init__()
      self.x, self.y, self.z, self.errorz = args
      self.ncount = 0
      
   def testfit(self):
      gMinuit = TMinuit(5)
      gMinuit.SetFCN(self.fcn)

      arglist = array( 'd', 10*[0.] )
      ierflg = ROOT.Long(1982)

      arglist[0] = 1
      gMinuit.mnexcm( "SET ERR", arglist, 1, ierflg )

    # Set starting values and step sizes for parameters
      vstart = array( 'd', ( 3,  1,  0.1,  0.01  ) )
      step   = array( 'd', ( 0.1, 0.1, 0.01, 0.001 ) )
      gMinuit.mnparm( 0, "a1", vstart[0], step[0], 0, 0, ierflg )
      gMinuit.mnparm( 1, "a2", vstart[1], step[1], 0, 0, ierflg )
      gMinuit.mnparm( 2, "a3", vstart[2], step[2], 0, 0, ierflg )
      gMinuit.mnparm( 3, "a4", vstart[3], step[3], 0, 0, ierflg )

    # Now ready for minimization step
      arglist[0] = 500
      arglist[1] = 1.
      gMinuit.mnexcm( "MIGRAD", arglist, 2, ierflg )

    # Print results
      # amin, edm, errdef = 0.18, 0.19, 0.20
      # nvpar, nparx, icstat = 1983, 1984, 1985
      # gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
      # gMinuit.mnprin( 3, amin )


   def fcn(self, npar, gin, f, par, iflag ):
      nbins = 5

    # calculate chisquare
      chisq, delta = 0., 0.
      for i in range(nbins):
         delta  = (self.z[i]-func(self.x[i], self.y[i],par)) / self.errorz[i]
         chisq += delta*delta

      f[0] = chisq
      self.ncount += 1


def main():
   z = array( 'f', ( 1., 0.96, 0.89, 0.85, 0.78 ) )
   errorz = array( 'f', 5*[0.01] )

   x = array( 'f', ( 1.5751, 1.5825,  1.6069,  1.6339,   1.6706  ) )
   y = array( 'f', ( 1.0642, 0.97685, 1.13168, 1.128654, 1.44016 ) )


   check = MinuitCheck(x, y, z, errorz)
   check.testfit()


if __name__ == '__main__':
   main()
