import os
import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    # settings
    inputfiles = sys.argv[1:]
    treename = 'tree'
    outputdir = 'output_plots'

    categories = {
      'recojet_isB': ('b-jets', 'darkviolet'),
      'recojet_isC': ('c-jets', 'mediumpurple'),
      'recojet_isS': ('s-jets', 'dodgerblue'),
      'recojet_isUD': ('ud-jets', 'deepskyblue'),
      #'recojet_isData': ('data', None)
    }

    variables = [
      'event_ntracks',
      'event_nselectedtracks',
      'event_nprimarytracks',
      'event_nsecondarytracks',
      'event_nsv',
      'event_nv0candidates',

      'recojet_pt',
      'recojet_eta',
      'recojet_theta',
      'recojet_phi',
      'recojet_e',
      'recojet_mass',
      'recojet_ntracks',
      'recojet_nselectedtracks',
      'recojet_nsecondarytracks',
      
      'nconst',
      'nphotons',
      'nchargedhad',
      'nneutralhad',
      'nel',
      'nmu',
      'recojet_nsv',
      'recojet_nv0candidates',
      'recojet_nkscandidates',
      'recojet_nlambdacandidates',

      #'sv_xrel',
      #'sv_yrel',
      #'sv_zrel',
      #'sv_thetarel',
      #'sv_phirel',
      #'sv_p',
      #'sv_prel',
      #'sv_chi2',
      #'sv_chi2Normalized',
      #'sv_ndof',
      #'sv_nTracks',
      #'sv_mass',
      #'sv_dxy',
      #'sv_dxyz',
      #'sv_cosPointing',
      #'sv_correctedMass',

      'v0cand_pdgId',
      'v0cand_xrel',
      'v0cand_yrel',
      'v0cand_zrel',
      'v0cand_thetarel',
      'v0cand_phirel',
      'v0cand_p',
      'v0cand_prel',
      'v0cand_chi2',
      'v0cand_chi2Normalized',
      'v0cand_ndof',
      'v0cand_nTracks',
      'v0cand_mass',
      'v0cand_dxy',
      'v0cand_dxyz',
      'v0cand_cosPointing',
      'v0cand_correctedMass',

      'pfcand_pt',
      #'pfcand_e',
      #'pfcand_ptrel_log',
      #'pfcand_erel_log',
      #'pfcand_thetarel',
      #'pfcand_phirel',
      #'pfcand_charge',

      #'pfcand_d0_wrt0',
      #'pfcand_z0_wrt0',
      #'pfcand_dxy',
      #'pfcand_dz',

      #'pfcand_Sip2dVal',
      #'pfcand_Sip2dSig',
      #'pfcand_Sip3dVal',
      #'pfcand_Sip3dSig',
      #'pfcand_JetDistVal',
      #'pfcand_JetDistSig',

      #'pfcand_linearSignedIP3D',
      #'pfcand_linearSignedIP3DSig',
      #'pfcand_transverseJetDistance',
      #'pfcand_longitudinalJetDistance',
        
      #'pfcand_dxydxy',
      #'pfcand_dzdz',

      #'pfcand_dEdx_pads_type',
      #'pfcand_dEdx_pads_value',
      #'pfcand_dEdx_pads_error',
      #'pfcand_dEdx_wires_type',
      #'pfcand_dEdx_wires_value',
      #'pfcand_dEdx_wires_error',      
    ]

    # make output dir if needed
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # read input files
    batches = []
    branches_to_read = list(categories.keys()) + variables
    for idx, inputfile in enumerate(inputfiles):
        print(f'Reading file {idx+1} / {len(inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches_to_read))
    events = ak.concatenate(batches)
    print(f'Read {len(events)} entries.')

    # make kinematic mask (optional)
    #kinematic_mask = np.ones(len(events)).astype(bool)
    kinematic_mask = (
        (events['recojet_pt'].to_numpy()>10)
        & (np.abs(events['recojet_eta'].to_numpy())<0.9)
        & (events['nconst'].to_numpy()>5)
    ).astype(bool)
    print(f'Made kinematic mask with {np.sum(kinematic_mask)} / {len(kinematic_mask)} entries passing.')

    # make category masks
    category_masks = {}
    for category in categories.keys():
        category_masks[category] = events[category].to_numpy().astype(bool)
    print(f'Found following number of entries per category:')
    for category in categories.keys(): print(f'  - {category}: {np.sum(category_masks[category])}')
    print(f'  -> total: {sum([np.sum(v) for v in category_masks.values()])}')

    # loop over variables to plot
    for variable in variables:
        print(f'Running on variable {variable}...')

        # get data
        data = {}
        for category in categories.keys():
            this_data = events[variable][(kinematic_mask) & (category_masks[category])]
            
            # strategies for flattening per-constituent data
            if variable.startswith('pfcand_'):

                strategy = 'leading' # choose from "flatten" or "leading"

                # approach 1: take all constituents
                if strategy=='flatten': this_data = ak.flatten(this_data)

                # approach 2: take leading constituent
                # note: constituents do not seem to be pt-ordered by default!
                elif strategy=='leading':
                    if len(this_data)==0: this_data = ak.Array([])
                    else:
                        len_mask = np.nonzero(ak.num(this_data))[0].to_numpy().astype(bool)
                        this_data = this_data[len_mask]
                        pt = events['pfcand_pt'][(kinematic_mask) & (category_masks[category])]
                        pt = pt[len_mask]
                        ids = ak.argmax(pt, axis=1).to_numpy()
                        ids = np.array(list(ids))
                        this_data = this_data[np.arange(len(this_data)), ids]

                else: raise Exception(f'Strategy {strategy} not recognized.')

            # strategries for flattening secondary vertex data
            if variable.startswith('sv_'):
                
                # approach 1: take all constituents
                this_data = ak.flatten(this_data)

            # strategies for flattening v0 candidate data
            if variable.startswith('v0cand_'):

                # approach 1: take all constituents
                this_data = ak.flatten(this_data)

            # parsing
            this_data = this_data.to_numpy()
            if np.isnan(this_data).any():
                msg = 'WARNING: replacing NaN by 0...'
                print(msg)
                np.nan_to_num(this_data, copy=False, nan=0)
            data[category] = this_data

        # optional: ignore dummy values
        # note: to clean up in upstream code; some dummies are set to -9 and others to -1...
        for category in categories.keys():
            this_data = data[category]
            mask = (np.abs(this_data+9)>1e-12).astype(bool)
            mask = ((mask) & (np.abs(this_data+1)>1e-12).astype(bool))
            data[category] = this_data[mask]

        # special case: ignore failed measurements for dEdx
        if variable == 'pfcand_dEdx_wires_value' or variable == 'pfcand_dEdx_pads_value':
            for category in categories.keys():
                this_data = data[category]
                mask = (this_data > 0.1).astype(bool)
                data[category] = this_data[mask]

        # group categories in single array
        data_array = np.concatenate(list(data.values()))

        # printouts for testing
        print(f'--- test output for variable {variable} ---')
        print(data_array)
        print('min: ', np.amin(data_array))
        print('max: ', np.amax(data_array))
        print('mean: ', np.mean(data_array))
        print('std: ', np.std(data_array))
        print('-----')

        # determine whether variable is integer
        is_integer = False
        if np.all(data_array - data_array.astype(int) < 1e-6): is_integer = True

        # determine suitable binning
        mask = (np.abs(data_array+9)>1e-3).astype(bool)
        npass = np.sum(mask.astype(int))
        if npass==0:
            msg = 'WARNING: no instances pass mask; skipping this variable.'
            print(msg)
            continue
        minv = np.quantile(data_array[mask], 0.05)
        maxv = np.quantile(data_array[mask], 0.95)
        if minv < 0:
            maxv = max(maxv, abs(minv))
            minv = -maxv
        if is_integer:
            maxv = max(maxv, 5)
            #minv = int(minv) - 0.5
            minv = -0.5
            maxv = int(maxv) + 0.5
        # special cases (hard-coded)
        if variable == 'sv_chi2Normalized': (minv, maxv) = (-2, 25)
        num_edges = 51
        if is_integer: num_edges = int(maxv - minv) + 1
        bins = np.linspace(minv, maxv, num=num_edges)

        # make figure
        fig, ax = plt.subplots()
        for category, settings in categories.items():
            hist = np.histogram(data[category], bins=bins)[0]
            errors = np.sqrt(hist)
            binwidths = bins[1:] - bins[:-1]
            integral = np.sum(np.multiply(hist, binwidths))
            if integral > 0:
                hist = hist / integral
                errors = errors / integral
            ax.stairs(hist, edges=bins, linewidth=2, label=settings[0], color=settings[1])
            ax.stairs(hist+errors, baseline=hist-errors, edges=bins, fill=True, color=settings[1], alpha=0.2)

        # plot aesthetics
        ax.set_ylabel('Jets (normalized)')
        ax.set_xlabel(variable)
        ax.legend()
        ax.grid(which='both', axis='both')

        # save figure
        fig.tight_layout()
        outputfile = os.path.join(outputdir, variable+'.png')
        fig.savefig(outputfile)

        # same with log scale
        ax.set_yscale('log')
        fig.tight_layout()
        outputfile = os.path.join(outputdir, variable+'_log.png')
        fig.savefig(outputfile)

        # close figures to save memory
        plt.close()
