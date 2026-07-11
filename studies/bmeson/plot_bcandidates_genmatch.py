# Check the semileptonic B -> D* l nu tag candidates (BCandidates_*), split by
# gen-match and by lepton charge sign ("right-sign" vs "wrong-sign").
#
# This requires stage-1 ntuples produced with do_bmeson_candidates = True on
# simulation (BCandidates_hasGenMatch is a dummy False for data, same as for
# DStarCandidates_hasGenMatch).
#
# Purpose: this is primarily meant to empirically check the assumed "right-sign"
# lepton charge convention (lepton charge opposite the D* candidate's kaon charge),
# which was recalled from memory rather than pulled from a live reference when
# implementing analyzers/analyzer_bmesonfinder.cxx -- if the convention is correct,
# gen-matched candidates should be predominantly right-sign. Secondarily, it's a
# first look at the visible (D* + lepton) mass distribution and overall purity, same
# spirit as studies/dstar/plot_dstar_genmatch.py for the D* candidates themselves.
#
# Note on the gen-matching itself: same "poor man's" geometric match as for the D*
# candidates (see analyzer_dstarfinder.cxx / analyzer_bmesonfinder.cxx docstrings),
# matched against the D* candidate's own direction (not the visible D*+lepton system).


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
    outputfile = 'bcandidates_genmatch_signcheck.png'
    bins = np.linspace(0, 6, 41)

    if len(inputfiles)==0:
        raise Exception('Please provide at least one input file (stage-1 ntuple).')

    # read input files
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}...')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(['BCandidates_mass', 'BCandidates_isRightSign', 'BCandidates_hasGenMatch']))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} events.')

    # flatten to one entry per B tag candidate
    mass = ak.to_numpy(ak.flatten(events['BCandidates_mass']))
    rightsign = ak.to_numpy(ak.flatten(events['BCandidates_isRightSign']))
    matched = ak.to_numpy(ak.flatten(events['BCandidates_hasGenMatch']))
    print(f'Found {len(mass)} B tag candidates'
          + f' ({matched.sum()} gen-matched, {(~matched).sum()} not matched).')

    # the actual sign-convention check: among gen-matched candidates, what fraction
    # are right-sign? should be close to 1 if the assumed convention is correct.
    if matched.sum() > 0:
        print(f'Among gen-matched candidates: {rightsign[matched].sum()} right-sign,'
              + f' {(~rightsign[matched]).sum()} wrong-sign'
              + f' ({100*rightsign[matched].mean():.1f}% right-sign).')
    if (~matched).sum() > 0:
        print(f'Among non-gen-matched candidates: {rightsign[~matched].sum()} right-sign,'
              + f' {(~rightsign[~matched]).sum()} wrong-sign'
              + f' ({100*rightsign[~matched].mean():.1f}% right-sign).')

    # make one panel for right-sign, one for wrong-sign, each showing the
    # matched/non-matched split of the visible (D* + lepton) mass
    categories = {
        'right-sign': rightsign,
        'wrong-sign': ~rightsign,
    }

    fig, axs = plt.subplots(1, len(categories), figsize=(6*len(categories), 5), sharey=True)
    for ax, (label, mask) in zip(axs, categories.items()):
        cand_mass = mass[mask]
        cand_matched = matched[mask]
        ax.hist(
            [cand_mass[~cand_matched], cand_mass[cand_matched]],
            bins=bins, stacked=True,
            color=['lightgray', 'dodgerblue'],
            label=['not gen-matched', 'gen-matched (proxy)'],
            edgecolor='black', linewidth=0.5,
        )
        ax.set_xlabel('Visible (D* + lepton) mass [GeV]', fontsize=13)
        ax.legend(title=label, fontsize=12,
                  title_fontproperties={'weight': 'bold', 'size': 13})
    axs[0].set_ylabel('Candidates / bin', fontsize=13)

    fig.tight_layout()
    fig.savefig(outputfile, dpi=150)
    print(f'Output written to {outputfile}.')
