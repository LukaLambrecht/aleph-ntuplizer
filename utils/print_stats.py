# Print some basic data statistics such as:
# - number of events
# - number of runs
# - number of fills

# Note: assumed to run on the raw input files,
# before processing!


import os
import sys
import uproot
import numpy as np
import awkward as ak


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'
    branches_to_read = [
        'EventHeader.eventNumber',
        'EventHeader.runNumber',
        #'EventHeader.fillNumber'
    ]

    # read input files
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches_to_read))
    events = ak.concatenate(batches)

    # get number of events
    nevents = len(events)

    # get number of runs
    runs = np.unique(events['EventHeader.runNumber'])
    nruns = len(runs)

    # get number of events per run
    nevents_per_run = {}
    for run in runs:
        nevents_per_run[run] = np.sum((events['EventHeader.runNumber']==run))
    nevents_per_run_avg = np.mean(list(nevents_per_run.values()))
    nevents_per_run_std = np.std(list(nevents_per_run.values()))

    # print some info
    print(f'Number of events: {nevents}')
    print(f'Number of runs: {nruns}')
    print(f'Number of events per run: {nevents_per_run_avg} +- {nevents_per_run_std}')
