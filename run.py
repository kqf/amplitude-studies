#!/usr/bin/python

from amplitude.amplitude import DataManager

def main():
    m = DataManager('alldata_v1_4.dat')
    m.plot_approximation()
    raw_input('Press any key ...')

if __name__ == "__main__":
    main()
