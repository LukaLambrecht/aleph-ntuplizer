# Plot D0 candidate properties (mass, vertex distance to the primary vertex), split
# into gen-matched and non-matched candidates, with separate panels per gen-level
# event flavour (bb, cc, other/light).
#
# This requires stage-1 ntuples produced with do_dzero_candidates = True on simulation
# (DZeroCandidates_hasGenMatch is a dummy False for data, so running this on data
# just puts everything in the "non-matched" category; genEventType is also a dummy
# -1 for data, so it will show up entirely in the "other" category below).
#
# Note on the gen-matching itself: it is a "poor man's" geometric match (is there a
# gen-level D0 within dR of the reconstructed D0 candidate?), not a decay-chain match,
# because mother-daughter links are not populated in this dataset's gen-particle
# collection (checked directly: MCParticles parents/daughters begin/end are zero-length
# for every particle). So this cannot be used to tune selection cuts against individual
# decay products, and a "matched" candidate could in principle be an accidental nearby
# gen D0 rather than the true parent -- but it's still a reasonable first-order check
# of how much of the reconstructed peak is real vs combinatorial fakes.
#
# Note on the flavour categories: genEventType is the (absolute) PDG ID of the first
# quark found in the gen-particle collection (see GenEventType::get_genEventTypeFromFirstQuark
# in analyzers/analyzer_geneventtype.cxx), i.e. 1-3 for d/u/s (light), 4 for c, 5 for b.
# "other (light)" below groups 1-3 together, plus anything not in {4, 5} as a catch-all
# (e.g. failed classification, or data), so the three panels always add up to the total.
#
# Note on the same-sign overlay: DZeroCandidates_isOppositeSign is False for K/pi
# candidates built from a same-charge pair (kept purely as a combinatorial background
# estimate, see analyzer_dzerofinder.cxx). It is drawn as an unfilled step histogram on
# the same axes as the (opposite-sign) signal candidates, not split by gen-match
# (gen-matching a same-charge fake is not meaningful) and not normalized to anything --
# just the raw same-sign yield, for a direct visual comparison of shape and level
# against the opposite-sign peak.
#
# Note on the D0 vertex distance: DZeroCandidates_vtx_dxyz is the 3D distance between
# the fitted D0 (K, pi) vertex and the primary vertex, i.e. how far the D0 travelled
# before decaying. This is used to separate almost-prompt D0 decays (from c-jets) from
# displaced D0 decays (where the D0 is a decay product of a B meson that itself
# travelled some distance first).
#
# Note on units: despite various code comments in analyzer_svfinder.cxx / the FCCAnalyses
# VertexingUtils helpers claiming vertex positions are "in mm", they are actually in cm for
# this dataset. analyzer_tracktools.cxx has a (deliberately disabled) cm->mm conversion for
# the raw D0/Z0 track parameters ("keep everything in cm here... to stay in sync with earlier
# version"), and the vertex fit (VertexFitterMod) is told Units_mm=true, which in FCCAnalyses
# just means "don't rescale the input" -- so cm in, cm out. Cross-checked empirically: the
# reconstructed beam-spot spread (GenPV_x/y std) only matches the known ALEPH/LEP beam size
# (sigma_x ~ 150 um, sigma_y ~ 6 um) if the raw numbers are interpreted as cm, not mm.
#
# Note on why this mass peak looks much less clean than the D* - D0 mass difference
# peak in plot_dstar.py: see the discussion in that conversation -- in short, the D*
# mass difference benefits from a large cancellation of correlated track-momentum
# measurement errors shared between the D0 and D* masses (both are built from the same
# K/pi tracks), which is not available here since there is no companion soft pion to
# subtract against. The D* requirement was also acting as a strong background-rejection
# cut in its own right (very few random K/pi combinatorics also have a correlated extra
# track at the right D*-D0 Q-value), which is likewise absent for the bare D0 search.


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
    branches = ['genEventType', 'DZeroCandidates_mass', 'DZeroCandidates_hasGenMatch',
                'DZeroCandidates_isOppositeSign', 'DZeroCandidates_vtx_dxyz']
    batches = []
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}...')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} events.')

    # broadcast the per-event genEventType to the per-candidate (jagged) structure,
    # then flatten everything to one entry per D0 candidate
    genEventType_bcast, _ = ak.broadcast_arrays(events['genEventType'], events['DZeroCandidates_mass'])
    mass = ak.to_numpy(ak.flatten(events['DZeroCandidates_mass']))
    vtx_dxyz = ak.to_numpy(ak.flatten(events['DZeroCandidates_vtx_dxyz']))
    matched = ak.to_numpy(ak.flatten(events['DZeroCandidates_hasGenMatch']))
    oppositesign = ak.to_numpy(ak.flatten(events['DZeroCandidates_isOppositeSign']))
    evttype = ak.to_numpy(ak.flatten(genEventType_bcast))
    print(f'Found {len(mass)} D0 candidates'
          + f' ({oppositesign.sum()} opposite-sign, {(~oppositesign).sum()} same-sign);'
          + f' among opposite-sign: {matched[oppositesign].sum()} gen-matched,'
          + f' {(~matched[oppositesign]).sum()} not matched.')

    # make mass plot
    print('Making mass plot...')
    make_stacked_plot(
        mass, matched, oppositesign, evttype,
        bins=np.linspace(1.7, 2.0, 61),
        xlabel='D0 mass [GeV]',
        outputfile='dzero_mass.png',
        refline=(1.86484, 'PDG D0 mass'),
    )

    # make D0 vertex - primary vertex distance plot
    print('Making D0 vertex distance plot...')
    make_stacked_plot(
        vtx_dxyz, matched, oppositesign, evttype,
        bins=np.linspace(0, 1.0, 51),
        xlabel='D0 vertex - primary vertex distance [cm]',
        outputfile='dzero_vtxdistance.png',
    )
