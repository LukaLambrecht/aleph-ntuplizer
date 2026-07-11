# Plot the D* - D0 mass difference, split into gen-matched and non-matched candidates,
# with separate panels per gen-level event flavour (bb, cc, other/light).
#
# This requires stage-1 ntuples produced with do_dstar_candidates = True on simulation
# (DStarCandidates_hasGenMatch is a dummy False for data, so running this on data
# just puts everything in the "non-matched" category; genEventType is also a dummy
# -1 for data, so it will show up entirely in the "other" category below).
#
# Note on the gen-matching itself: it is a "poor man's" geometric match (is there a
# gen-level D*(2010) within dR of the reconstructed D* candidate?), not a decay-chain
# match, because mother-daughter links are not populated in this dataset's gen-particle
# collection (checked directly: MCParticles parents/daughters begin/end are zero-length
# for every particle). So this cannot be used to tune selection cuts against individual
# decay products, and a "matched" candidate could in principle be an accidental nearby
# gen D* rather than the true parent -- but it's still a reasonable first-order check
# of how much of the reconstructed peak is real vs combinatorial fakes.
#
# Note on the flavour categories: genEventType is the (absolute) PDG ID of the first
# quark found in the gen-particle collection (see GenEventType::get_genEventTypeFromFirstQuark
# in analyzers/analyzer_geneventtype.cxx), i.e. 1-3 for d/u/s (light), 4 for c, 5 for b.
# "other (light)" below groups 1-3 together, plus anything not in {4, 5} as a catch-all
# (e.g. failed classification, or data), so the three panels always add up to the total.


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
    outputfile = 'dstar_genmatch_massdiff.png'
    bins = np.linspace(0.14, 0.20, 61)

    if len(inputfiles)==0:
        raise Exception('Please provide at least one input file (stage-1 ntuple).')

    # read input files
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}...')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(['genEventType', 'DStarCandidates_d0_massDiff', 'DStarCandidates_hasGenMatch']))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} events.')

    # broadcast the per-event genEventType to the per-candidate (jagged) structure,
    # then flatten everything to one entry per D* candidate
    genEventType_bcast, _ = ak.broadcast_arrays(events['genEventType'], events['DStarCandidates_d0_massDiff'])
    massdiff = ak.to_numpy(ak.flatten(events['DStarCandidates_d0_massDiff']))
    matched = ak.to_numpy(ak.flatten(events['DStarCandidates_hasGenMatch']))
    evttype = ak.to_numpy(ak.flatten(genEventType_bcast))
    print(f'Found {len(massdiff)} D* candidates'
          + f' ({matched.sum()} gen-matched, {(~matched).sum()} not matched).')

    # define flavour categories (see note above on genEventType)
    categories = {
        'bb events': (evttype == 5),
        'cc events': (evttype == 4),
        'uu / dd / ss events': ~np.isin(evttype, [4, 5]),
    }

    # make one panel per category
    fig, axs = plt.subplots(1, len(categories), figsize=(6*len(categories), 5), sharey=True)
    for ax, (label, mask) in zip(axs, categories.items()):
        cand_massdiff = massdiff[mask]
        cand_matched = matched[mask]
        print(f'  {label}: {mask.sum()} candidates'
              + f' ({cand_matched.sum()} gen-matched, {(~cand_matched).sum()} not matched).')

        ax.hist(
            [cand_massdiff[~cand_matched], cand_massdiff[cand_matched]],
            bins=bins, stacked=True,
            color=['lightgray', 'dodgerblue'],
            label=['not gen-matched', 'gen-matched (proxy)'],
            edgecolor='black', linewidth=0.5,
        )
        ax.axvline(0.14543, color='red', linestyle='--', linewidth=1, label='PDG D*-D0 mass diff')
        ax.set_xlabel('D* - D0 mass difference [GeV]', fontsize=13)
        ax.legend(title=label, fontsize=12,
                  title_fontproperties={'weight': 'bold', 'size': 13})
    axs[0].set_ylabel('Candidates / bin', fontsize=13)

    fig.tight_layout()
    fig.savefig(outputfile, dpi=150)
    print(f'Output written to {outputfile}.')
