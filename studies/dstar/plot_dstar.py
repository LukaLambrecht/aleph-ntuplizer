# Plot D* candidate properties (mass difference, vertex distance to the primary
# vertex), split into gen-matched and non-matched candidates, with separate panels
# per gen-level event flavour (bb, cc, other/light).
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
#
# Note on the same-sign overlay: DStarCandidates_isOppositeSign is False for K/pi2
# candidates built from a same-charge pair (kept purely as a combinatorial background
# estimate, see analyzer_dstarfinder.cxx -- same idea as the right-sign/wrong-sign
# check for the B candidates, but for the D0 K-pi combination itself). It is drawn
# as an unfilled step histogram on the same axes as the (opposite-sign) signal
# candidates, not split by gen-match (gen-matching a same-charge fake is not
# meaningful) and not normalized to anything -- just the raw same-sign yield, for a
# direct visual comparison of shape and level against the opposite-sign peak.
#
# Note on the D* vertex distance: DStarCandidates_dstarvtx_dxyz is the 3D distance
# between the fitted D* (K, pi2, slow pi) vertex and the primary vertex, i.e. how far
# the D* travelled before decaying. This is used to separate almost-prompt D* decays
# (from c-jets, where the D* is produced directly at the primary vertex) from
# displaced D* decays (where the D* is a decay product of a B meson that itself
# travelled some distance first).


import os
import sys
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt


def make_stacked_plot(values, matched, oppositesign, evttype, bins,
                       xlabel, outputfile, refline=None):
    '''
    Make one panel per flavour category, each showing a stacked gen-matched /
    non-matched histogram of values (for opposite-sign candidates only), plus an
    unfilled step histogram of the same-sign background estimate.
    '''
    categories = {
        'bb events': (evttype == 5),
        'cc events': (evttype == 4),
        'uu / dd / ss events': ~np.isin(evttype, [4, 5]),
    }

    fig, axs = plt.subplots(1, len(categories), figsize=(6*len(categories), 5), sharey=True)
    for ax, (label, mask) in zip(axs, categories.items()):
        os_mask = mask & oppositesign
        ss_mask = mask & ~oppositesign
        cand_values = values[os_mask]
        cand_matched = matched[os_mask]
        print(f'  {label}: {os_mask.sum()} opposite-sign candidates'
              + f' ({cand_matched.sum()} gen-matched, {(~cand_matched).sum()} not matched),'
              + f' {ss_mask.sum()} same-sign candidates.')

        ax.hist(
            [cand_values[~cand_matched], cand_values[cand_matched]],
            bins=bins, stacked=True,
            color=['lightgray', 'dodgerblue'],
            label=['not gen-matched', 'gen-matched (proxy)'],
            edgecolor='black', linewidth=0.5,
        )
        ax.hist(
            values[ss_mask], bins=bins,
            histtype='step', color='black', linewidth=1.5,
            label='same-sign (background estimate)',
        )
        if refline is not None:
            ax.axvline(refline[0], color='red', linestyle='--', linewidth=1, label=refline[1])
        ax.set_xlabel(xlabel, fontsize=13)
        ax.legend(title=label, fontsize=12,
                  title_fontproperties={'weight': 'bold', 'size': 13})
    axs[0].set_ylabel('Candidates / bin', fontsize=13)

    fig.tight_layout()
    fig.savefig(outputfile, dpi=150)
    print(f'Output written to {outputfile}.')


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'events'

    if len(inputfiles)==0:
        raise Exception('Please provide at least one input file (stage-1 ntuple).')

    # read input files
    branches = ['genEventType', 'DStarCandidates_d0_massDiff', 'DStarCandidates_hasGenMatch',
                'DStarCandidates_isOppositeSign', 'DStarCandidates_dstarvtx_dxyz']
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}...')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} events.')

    # broadcast the per-event genEventType to the per-candidate (jagged) structure,
    # then flatten everything to one entry per D* candidate
    genEventType_bcast, _ = ak.broadcast_arrays(events['genEventType'], events['DStarCandidates_d0_massDiff'])
    massdiff = ak.to_numpy(ak.flatten(events['DStarCandidates_d0_massDiff']))
    dstarvtx_dxyz = ak.to_numpy(ak.flatten(events['DStarCandidates_dstarvtx_dxyz']))
    matched = ak.to_numpy(ak.flatten(events['DStarCandidates_hasGenMatch']))
    oppositesign = ak.to_numpy(ak.flatten(events['DStarCandidates_isOppositeSign']))
    evttype = ak.to_numpy(ak.flatten(genEventType_bcast))
    print(f'Found {len(massdiff)} D* candidates'
          + f' ({oppositesign.sum()} opposite-sign, {(~oppositesign).sum()} same-sign);'
          + f' among opposite-sign: {matched[oppositesign].sum()} gen-matched,'
          + f' {(~matched[oppositesign]).sum()} not matched.')

    # make mass difference plot
    print('Making mass difference plot...')
    make_stacked_plot(
        massdiff, matched, oppositesign, evttype,
        bins=np.linspace(0.14, 0.20, 61),
        xlabel='D* - D0 mass difference [GeV]',
        outputfile='dstar_massdiff.png',
        refline=(0.14543, 'PDG D*-D0 mass diff'),
    )

    # make D* vertex - primary vertex distance plot
    print('Making D* vertex distance plot...')
    make_stacked_plot(
        dstarvtx_dxyz, matched, oppositesign, evttype,
        bins=np.linspace(0, 1.0, 51),
        xlabel='D* vertex - primary vertex distance [mm]',
        outputfile='dstar_vtxdistance.png',
    )
