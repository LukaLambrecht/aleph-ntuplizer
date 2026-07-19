import os
import sys
import ROOT
import numpy as np
import awkward as ak
import uproot
import argparse


if __name__=='__main__':

    # read command line args
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', required=True)
    parser.add_argument('-o', '--outputfile', required=True)
    parser.add_argument('--entry_start', default=None)
    parser.add_argument('--entry_stop', default=None)
    args = parser.parse_args()

    # read input file
    print(f'Reading file {args.inputfile}...')
    df = ROOT.RDataFrame("events", args.inputfile)
    cols = [str(c) for c in df.GetColumnNames()]
    nentries = len(ak.from_rdataframe(df, columns=['eventNumber'])['eventNumber'])
    print(f'Read file {args.inputfile} with {nentries} entries {len(cols)} branches.')

    # define helper function for parsing number of entries to process
    def parse_nentries(nentries, max_nentries):
        if nentries is None: return 0
        nentries = float(nentries)
        if nentries.is_integer(): return int(nentries)
        return min(max_nentries, int(nentries * max_nentries))

    # limit number of entries to process
    if args.entry_start is not None or args.entry_stop is not None:
        entry_start = parse_nentries(args.entry_start, nentries)
        entry_stop = parse_nentries(args.entry_stop, nentries)
        if entry_stop < entry_start: entry_stop = nentries
        print(f'Will process entries {entry_start} - {entry_stop} (out of {nentries} available).')
        df = df.Range(entry_start, entry_stop)

    # filter events
    df = df.Filter("Event_njets >= 2")

    # loop over columns
    dout = {}
    for col in cols:
   
        # special case for some per-event variables that are written as length-1 vectors per event
        # (skip for now)
        if col in ['eventNumber', 'runNumber']: continue

        # special cases: skip some variables that were not written to the per-jet ntuples in the previous version
        # (only kept to be able to exactly sychronize with previous version; doesn't hurt to just add them in later)
        skipcols = [
            'Event_de', 'Event_dphi', 'Event_dpt', 'Event_dtheta', 'Event_mass',
            'Event_nConstituentsSum', 'Event_nRecoParticles', 'Event_nTracksPerJetSum',
            'Event_nTracksWithReco', 'Event_nTracksWithoutReco', 'Event_njets',
            'JetsConstituents_C', 'JetsConstituents_ct',
            'JetsConstituents_jetAxisPCAToTrack_x', 'JetsConstituents_jetAxisPCAToTrack_y', 'JetsConstituents_jetAxisPCAToTrack_z',
            'JetsConstituents_linePCAToPrimaryVertex_x',
            'JetsConstituents_linePCAToPrimaryVertex_y',
            'JetsConstituents_linePCAToPrimaryVertex_z',
            'JetsConstituents_trackPCAToJetAxis_x', 'JetsConstituents_trackPCAToJetAxis_y', 'JetsConstituents_trackPCAToJetAxis_z',
            'JetsConstituents_pdgId',
            'JetsConstituents_omega_wrt0', 'JetsConstituents_tanlambda_wrt0',
            'Jets_px', 'Jets_py', 'Jets_pz',

            # DZeroCandidates_* are a flat, event-level collection (one entry per
            # D0 candidate, not per jet); same reasoning as DStarCandidates_* below.
            'Event_nDZeroCandidates',
            'DZeroCandidates_mass', 'DZeroCandidates_pt', 'DZeroCandidates_eta', 'DZeroCandidates_phi',
            'DZeroCandidates_k_pt', 'DZeroCandidates_k_eta', 'DZeroCandidates_k_phi', 'DZeroCandidates_k_charge',
            'DZeroCandidates_pi_pt', 'DZeroCandidates_pi_eta', 'DZeroCandidates_pi_phi', 'DZeroCandidates_pi_charge',
            'DZeroCandidates_tr1tr2_deltaR', 'DZeroCandidates_vtx_chi2Normalized',
            'DZeroCandidates_vtx_dxy', 'DZeroCandidates_vtx_dxyz',
            'DZeroCandidates_isOppositeSign', 'DZeroCandidates_hasGenMatch',

            # DStarCandidates_* are a flat, event-level collection (one entry per
            # D* candidate, not per jet); they don't fit the "truncate/repeat to
            # the leading 2 jets" reshaping below, so skip them here rather than
            # have them silently mis-truncated/mis-zipped as if jet-indexed.
            'Event_nDStarCandidates',
            'DStarCandidates_mass', 'DStarCandidates_pt', 'DStarCandidates_eta', 'DStarCandidates_phi',
            'DStarCandidates_d0_mass', 'DStarCandidates_d0_pt', 'DStarCandidates_d0_eta', 'DStarCandidates_d0_phi',
            'DStarCandidates_d0_massDiff',
            'DStarCandidates_pi1_pt', 'DStarCandidates_pi1_eta', 'DStarCandidates_pi1_phi', 'DStarCandidates_pi1_charge',
            'DStarCandidates_k_pt', 'DStarCandidates_k_eta', 'DStarCandidates_k_phi', 'DStarCandidates_k_charge',
            'DStarCandidates_pi2_pt', 'DStarCandidates_pi2_eta', 'DStarCandidates_pi2_phi', 'DStarCandidates_pi2_charge',
            'DStarCandidates_tr1tr2_deltaR', 'DStarCandidates_tr3d0_deltaR',
            'DStarCandidates_d0vtx_chi2Normalized', 'DStarCandidates_dstarvtx_chi2Normalized',
            'DStarCandidates_d0vtx_dxy', 'DStarCandidates_d0vtx_dxyz',
            'DStarCandidates_dstarvtx_dxy', 'DStarCandidates_dstarvtx_dxyz',
            'DStarCandidates_isOppositeSign', 'DStarCandidates_hasGenMatch',

            # BCandidates_* are likewise a flat, event-level collection (one entry per
            # semileptonic B tag candidate), same reasoning as DStarCandidates_* above.
            'Event_nBCandidates',
            'BCandidates_dstar_idx',
            'BCandidates_lepton_pt', 'BCandidates_lepton_eta', 'BCandidates_lepton_phi',
            'BCandidates_lepton_charge', 'BCandidates_lepton_type',
            'BCandidates_mass', 'BCandidates_pt', 'BCandidates_eta', 'BCandidates_phi',
            'BCandidates_dstarlepton_deltaR', 'BCandidates_isRightSign',
            'BCandidates_bvtx_chi2Normalized', 'BCandidates_bvtx_dxy', 'BCandidates_bvtx_dxyz',
            'BCandidates_dstarHasGenMatch', 'BCandidates_hasGenMatch',
        ]
        if col in skipcols: continue
 
        # read values as awkward array
        val = ak.from_rdataframe(df, columns=[col])[col]
        depth = val.layout.minmax_depth[1]

        # limit to 2 jets per event
        if depth >= 2: val = val[:, :2]

        # expand per-event variables
        if depth == 1: val = np.repeat(val[:, np.newaxis], 2, axis=1)

        # parse key
        newkey = col
        newkey = newkey.replace("Event_", "event_")
        newkey = newkey.replace("Jets_", "recojet_")
        newkey = newkey.replace("JetsConstituents_", "pfcand_")
        newkey = newkey.replace("SecondaryVertices_", "sv_")
        newkey = newkey.replace("V0Candidates_", "v0cand_")
        newkey = newkey.replace("PV_", "pv_")

        # special cases for variable renaming
        # (only kept to be able to exactly synchronize with previous version
        # and avoid changing any downstream code)
        translation = {
        	'Event_Bz': 'bz',
        	'Event_nPrimaryTracks': 'event_nprimarytracks',
        	'Event_nSecondaryTracks': 'event_nsecondarytracks',
        	'Event_nSelectedTracks': 'event_nselectedtracks',
            'Event_nSV': 'event_nsv',
            'Event_nTracks': 'event_ntracks',
            'Event_nV0Candidates': 'event_nv0candidates',
        	'Jets_nChargedHad': 'nchargedhad',
        	'Jets_nNeutralHad': 'nneutralhad',
        	'Jets_nPhoton': 'nphotons',
        	'Jets_nEl': 'nel',
        	'Jets_nMu': 'nmu',
        	'Jets_nConstituents': 'nconst',
            'Jets_nV0Candidates': 'recojet_nv0candidates',
            'Jets_nKsCandidates': 'recojet_nkscandidates',
            'Jets_nLambdaCandidates': 'recojet_nlambdacandidates',
            'Jets_nSecondaryTracksPerJet': 'recojet_nsecondarytracks',
            'Jets_nSelectedTracksPerJet': 'recojet_nselectedtracks',
            'Jets_nTracksPerJet': 'recojet_ntracks',
            'Jets_nSV': 'recojet_nsv',
            'PV_x': 'pvx',
            'PV_y': 'pvy',
            'PV_z': 'pvz',
            'JetsConstituents_omega_tanlambda_cov': 'pfcand_cctgtheta',
            'JetsConstituents_omega_z0_cov': 'pfcand_cdz',
            'JetsConstituents_tanlambda_cov': 'pfcand_detadeta',
            'JetsConstituents_tanlambda_z0_cov': 'pfcand_dlambdadz',
            'JetsConstituents_phi0': 'pfcand_dphi',
            'JetsConstituents_phi0_cov': 'pfcand_dphidphi',
            'JetsConstituents_phi0_d0_cov': 'pfcand_dphidxy',
            'JetsConstituents_omega_cov': 'pfcand_dptdpt',
            'JetsConstituents_omega_d0_cov': 'pfcand_dxyc',
            'JetsConstituents_tanlambda_d0_cov': 'pfcand_dxyctgtheta',
            'JetsConstituents_d0_cov': 'pfcand_dxydxy',
            'JetsConstituents_d0_z0_cov': 'pfcand_dxydz',
            'JetsConstituents_z0_cov': 'pfcand_dzdz',
            'JetsConstituents_omega_phi0_cov': 'pfcand_phic',
            'JetsConstituents_tanlambda_phi0_cov': 'pfcand_phictgtheta',
            'JetsConstituents_phi0_z0_cov': 'pfcand_phidz'
            
        }
        for orig, new in translation.items():
            if col == orig: newkey = new

        # do flattening
        newval = ak.flatten(val, axis=1)

        # add to output dict
        dout[newkey] = newval

    # calculate jet flavour
    eventflavour = dout['genEventType'].to_numpy().astype(int)
    dout['recojet_flavour'] = eventflavour
    dout['recojet_isData'] = (eventflavour<0)
    dout['recojet_isQ'] = (eventflavour>0)
    dout['recojet_isT'] = (eventflavour==6)
    dout['recojet_isB'] = (eventflavour==5)
    dout['recojet_isC'] = (eventflavour==4)
    dout['recojet_isS'] = (eventflavour==3)
    dout['recojet_isU'] = (eventflavour==2)
    dout['recojet_isD'] = (eventflavour==1)
    dout['recojet_isG'] = (eventflavour==0)
    dout['recojet_isUD'] = ((eventflavour==1) | (eventflavour==2))
    dout['recojet_isUDG'] = ((eventflavour==0) | (eventflavour==1) | (eventflavour==2))
    dout['recojet_isUDS'] = ((eventflavour==1) | (eventflavour==2) | (eventflavour==3))
    dout['recojet_isUDSG'] = ((eventflavour>=0) & (eventflavour<=3))

    # zip to avoid counter branches to be created
    zipped = {}
    collections = ['pfcand_', 'sv_', 'v0cand_']
    for collection in collections:
        thisdict = {}
        collection_name = collection.strip('_')
        for key, val in dout.items():
            if not key.startswith(collection): continue
            newkey = key.replace(collection, '')
            thisdict[newkey] = val
        if len(thisdict) == 0: continue
        zipped[collection_name] = ak.zip(thisdict)
    copies = [
        'bz',
        'nchargedhad', 'nneutralhad', 'nphotons', 'nel', 'nmu', 'nconst',
        'pvx', 'pvy', 'pvz',
        'event', 'recojet'
    ]
    for key, val in dout.items():
        for copy in copies:
            if key.startswith(copy):
                zipped[key] = val
                break

    # write output file
    with uproot.recreate(args.outputfile) as fout:
        fout["tree"] = zipped
    print(f'Output file {args.outputfile} written.')
