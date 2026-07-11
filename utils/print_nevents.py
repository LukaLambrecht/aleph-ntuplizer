# Simple utility script to print the number of events in a root file


import os
import sys
import ROOT
import numpy as np


if __name__=='__main__':

    # read input files
    inputfiles = sys.argv[1:]

    # get number of events per file
    nevents= {}
    nfiles = len(inputfiles)
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {nfiles}...', end='\r')
        f = ROOT.TFile.Open(inputfile)
        try:
            # for raw input files in EDM4HEP format
            tree = f.Get("events")
            tree.GetEntries()
        except:
            # for ntuples
            tree = f.Get("tree")
        nevents[inputfile] = tree.GetEntries()

    # print results
    print('Found following number of events:')
    for key, val in nevents.items():
        print(f'  - {key}: {val}')
    print('Summary statistics:')
    ns = np.array(list(nevents.values()))
    ntot = np.sum(ns)
    navg = np.mean(ns)
    nstd = np.std(ns)
    nmax = np.amax(ns)
    nmin = np.amin(ns)
    print(f' - Total: {ntot}')
    print(f' - Average: {navg}')
    print(f' - Std: {nstd}')
    print(f' - Max: {nmax}')
    print(f' - Min: {nmin}')
