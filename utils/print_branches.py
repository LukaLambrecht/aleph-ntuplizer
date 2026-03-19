# Simple utility script to print the available branches in a root file


import os
import sys
import ROOT


if __name__=='__main__':

    inputfile = sys.argv[1]
    print(f'Reading input file {inputfile}...')
    f = ROOT.TFile.Open(inputfile)
    tree = f.Get("events")
    for branch in tree.GetListOfBranches():
        print(f' - {branch.GetName()}: {branch.GetClassName()}')
