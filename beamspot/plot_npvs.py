# Plot number of primary vertices per run

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    # get input file
    inputfile = sys.argv[1]

    # read input file and get number of primary vertices
    with open(inputfile, 'r') as f:
        data = json.load(f)
    nevents = np.array([val['nevents'] for val in data.values()])
    npvs = np.array([val['npvs'] for val in data.values()])

    # clip
    amax = 7000
    nevents = np.clip(nevents, a_min=None, a_max=amax)
    npvs = np.clip(npvs, a_min=None, a_max=amax)

    # plot distributions
    fig, ax = plt.subplots()
    bins = np.linspace(-0.5, amax+0.5, num=50)
    # nevents
    counts = np.histogram(nevents, bins=bins)[0]
    errors = np.sqrt(counts)
    ax.stairs(counts+errors, baseline=counts-errors, edges=bins,
       color='dodgerblue', fill=True, alpha=0.3)
    ax.stairs(counts, edges=bins,
      color='dodgerblue', linewidth=2,
      label='Number of events')
    # npvs
    counts = np.histogram(npvs, bins=bins)[0]
    errors = np.sqrt(counts)
    ax.stairs(counts+errors, baseline=counts-errors, edges=bins,
       color='darkviolet', fill=True, alpha=0.3)
    ax.stairs(counts, edges=bins,
      color='darkviolet', linewidth=2,
      label='Number of events with valid primary vertex')
    ax.legend(fontsize=12)
    ax.set_ylabel('Runs', fontsize=12)
    ax.set_xlabel('Number of events', fontsize=12)

    # print median
    print(f'Median: {np.median(npvs)}')
    print(f'Mean: {np.mean(npvs)}')

    # save figure
    fig.tight_layout()
    fig.savefig('test.png')
