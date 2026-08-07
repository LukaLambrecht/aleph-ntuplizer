# Plot the reconstructed mass of V0 candidates (Ks and Lambda), showing all
# candidates vs. tight-only candidates as separate lines.
#
# This requires stage-1 ntuples produced with do_v0_candidates = True (either
# finder). V0Candidates_tight is -1 (not applicable) for every candidate when
# do_new_v0_finder = False, so with the old finder the "tight only" line will
# be empty; this plot is primarily meant for the new two-tier finder
# (analyzer_v0new.cxx, do_new_v0_finder = True), where V0Candidates_tight is
# 1 for candidates passing the adopted tight package and 0 for the loose
# (ML-training) tier only.


import os
import sys
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'
    outputfile = 'v0_mass.png'
    # in older samples, the V0Candidates_tight branch is not available;
    # hence disable it here (maybe deduce automatically later).
    do_tight = False

    if len(inputfiles)==0:
        raise Exception('Please provide at least one input file (stage-1 ntuple).')

    # read input files
    branches = ['V0Candidates_pdgId', 'V0Candidates_mass']
    if do_tight: branches += ['V0Candidates_tight']
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}...')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} events.')

    # flatten to one entry per V0 candidate (V0Candidates_* is nested per jet,
    # i.e. event -> jet -> candidate, so flatten fully)
    pdgId = ak.to_numpy(ak.flatten(events['V0Candidates_pdgId'], axis=None))
    mass = ak.to_numpy(ak.flatten(events['V0Candidates_mass'], axis=None))
    tight = np.ones(mass.shape) * -1
    if do_tight: tight = ak.to_numpy(ak.flatten(events['V0Candidates_tight'], axis=None))
    print(f'Found {len(mass)} V0 candidates ({(tight==1).sum()} tight,'
          + f' {(tight==0).sum()} loose-only, {(tight==-1).sum()} n/a).')

    # one panel per hypothesis
    categories = {
        'Ks': (310, 0.497611, np.linspace(0.40, 0.60, 41)),
        'Lambda': (3122, 1.115683, np.linspace(1.08, 1.20, 41)),
    }

    fig, axs = plt.subplots(1, len(categories), figsize=(6*len(categories), 5))
    for ax, (label, (pdg, pdgmass, bins)) in zip(axs, categories.items()):
        sel = (pdgId == pdg)
        cand_mass = mass[sel]
        cand_tight = tight[sel] == 1
        print(f'  {label}: {sel.sum()} candidates ({cand_tight.sum()} tight)')

        ax.hist(
            cand_mass, bins=bins,
            histtype='step', color='dodgerblue', linewidth=1.5,
            label='all candidates',
        )
        ax.hist(
            cand_mass[cand_tight], bins=bins,
            histtype='step', color='darkorange', linewidth=1.5,
            label='tight only',
        )
        ax.axvline(pdgmass, color='red', linestyle='--', linewidth=1, label=f'PDG {label} mass')
        ax.set_xlabel(f'{label} candidate mass [GeV]', fontsize=13)
        ax.set_ylabel('Candidates / bin', fontsize=13)
        ax.legend(fontsize=12)

    fig.tight_layout()
    fig.savefig(outputfile, dpi=150)
    print(f'Output written to {outputfile}.')
