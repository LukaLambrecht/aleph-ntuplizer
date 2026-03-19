import os
import sys
import json
import uproot
import argparse
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fit_beamspot_simple(values):
    # do simple fit of mean and width of beamspot
    
    # remove outliers
    lower_bound = np.quantile(values, 0.01)
    upper_bound = np.quantile(values, 0.99)
    values = values[((values>lower_bound) & (values<upper_bound))]

    # calculate mean and std
    mean = np.mean(values)
    width = np.std(values)
    meanunc = width/np.sqrt(len(values))

    return (mean, width, meanunc)


def fit_beamspot_gaussian(values):
    # estimate mean and width by fitting a gaussian

    # remove outliers
    lower_bound = np.quantile(values, 0.01)
    upper_bound = np.quantile(values, 0.99)
    values = values[((values>lower_bound) & (values<upper_bound))]

    # binning
    bins = np.linspace(np.amin(values), np.amax(values), num=100)
    x = (bins[:-1] + bins[1:])/2
    y = np.histogram(values, bins=bins)[0]

    # do fit
    def gauss(x, a, mu, sigma):
        return a*np.exp(-0.5*np.square(np.divide((x-mu), sigma)))
    a0 = np.amax(y)
    mu0 = np.mean(values)
    sigma0 = np.std(values)
    parameters, cov = curve_fit(gauss, x, y, p0=(a0, mu0, sigma0))
    
    # get parameters
    mean = parameters[1]
    width = parameters[2]
    meanunc = np.sqrt(cov[1,1])

    return (mean, width, meanunc)


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfiles', required=True, nargs='+',
      help='Input files, usually the output of the previous step (run).')
    parser.add_argument('-o', '--outputdir', required=True,
      help='Output directory.')
    parser.add_argument('--do_run_plots', default=False, action='store_true',
      help='Make per-run plots of distributions of PV coordinates')
    parser.add_argument('--sim', default=False, action='store_true',
      help='Run in simulation mode with fake run numbers')
    args = parser.parse_args()

    # fixed settings
    treename = 'events'
    pv_vars = ['PV_x', 'PV_y', 'PV_z']
    id_vars = ['runNumber']

    # make output directory
    if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

    # read input files
    batches = []
    branches_to_read = pv_vars + id_vars
    for idx, inputfile in enumerate(args.inputfiles):
        print(f'Reading file {idx+1} / {len(args.inputfiles)}', end='\r')
        readstr = ':'.join([inputfile, treename])
        with uproot.open(readstr) as f:
            batches.append(f.arrays(branches_to_read))
    events = ak.concatenate(batches)

    # simulation mode: divide in fake runs
    if args.sim:
        run_size = 800 # typically about 800 events with valid PV per run in data
        run_numbers = []
        counter = 1
        while len(run_numbers)*run_size < len(events):
            run_numbers.append(np.ones(run_size)*counter)
            counter += 1
        run_numbers = np.concatenate(run_numbers).astype(int)
        run_numbers = run_numbers[:len(events)]
        events['runNumber'] = run_numbers

    # make mask for events with valid primary vertex
    pv_mask = (
      (np.abs(events['PV_x'])>1e-12)
      | (np.abs(events['PV_y'])>1e-12)
      | (np.abs(events['PV_z'])>1e-12)
    ).to_numpy().astype(bool)

    # find runs
    runs = np.unique(events['runNumber'].to_numpy())

    # loop over runs
    data = {}
    for runidx, run in enumerate(runs):
        data[run] = {}
        print(f'Processing run {run} ({runidx+1} / {len(runs)})...', end='\r')

        # make a mask for this run
        run_mask = np.squeeze((events['runNumber']==run).to_numpy()).astype(bool)
        tot_mask = (pv_mask & run_mask)

        # store number of events and number of events with valid primary vertex
        data[run]['nevents'] = int(np.sum(run_mask.astype(int)))
        data[run]['npvs'] = int(np.sum(tot_mask.astype(int)))

        # safety for too small runs
        if np.sum(tot_mask.astype(int)) < 100:
            msg = f'WARNING: skipping run {run} because too few events.'
            print(msg)
            data[run]['fits'] = None
            continue

        # loop over variables
        data[run]['fits'] = {}
        data_for_plotting = {}
        for varname in pv_vars:
            values = events[varname][tot_mask].to_numpy()
            lower_bound = np.quantile(values, 0.01)
            upper_bound = np.quantile(values, 0.99)
            values_for_plotting = values[((values>lower_bound) & (values<upper_bound))]
            data_for_plotting[varname] = values_for_plotting

            # calculate the central value and width
            #data[run][varname] = fit_beamspot_simple(values)
            data[run]['fits'][varname] = fit_beamspot_gaussian(values)

        if not args.do_run_plots: continue

        # make figure
        fig, axs = plt.subplots(nrows=3)
        for idx, varname in enumerate(pv_vars):
            ax = axs[idx]
            values = data_for_plotting[varname]
            mean = data[run][varname][0]
            width = data[run][varname][1]
            ax.hist(values, bins=50, histtype='step', linewidth=2)
            ax.axvline(x=mean, color='blue', linestyle='dashed')
            ax.axvline(x=mean+width, color='red', linestyle='dashed')
            ax.axvline(x=mean-width, color='red', linestyle='dashed')

            # plot aesthetics
            ax.set_ylabel('Events')
            ax.text(0.98, 0.95, varname, ha='right', va='top', transform=ax.transAxes)
            ax.text(0.98, 0.75, 'Mean: {:.2e} cm'.format(mean), ha='right', va='top', transform=ax.transAxes)
            ax.text(0.98, 0.55, 'Width: {:.2e} cm'.format(width), ha='right', va='top', transform=ax.transAxes)

        # more plot aesthetics
        axs[-1].set_xlabel('[cm]')
        
        axs[0].text(0.99, 1.05, f'Run {run}', ha='right', transform=axs[0].transAxes)

        # save figure
        fig.tight_layout()
        outputfile = os.path.join(args.outputdir, f'run_{run}.png')
        fig.savefig(outputfile)

        # close figures to save memory
        plt.close()
    print()

    # fill runs for which no sensible measurement could be made with sensible values
    for runidx, run in enumerate(runs):
        if data[run]['fits'] is None:
            if runidx==0:
                nextrunidx = runidx + 1
                nextrun = runs[nextrunidx]
                while data[nextrun]['fits'] is None:
                    nextrunidx += 1
                    nextrun = runs[nextrunidx]
                data[run]['fits'] = data[nextrun]['fits']
            else:
                previousrunidx = runidx - 1
                previousrun = runs[previousrunidx]
                while data[previousrun]['fits'] is None:
                    previousrunidx -= 1
                    previousrun = runs[previousrunidx]
                data[run]['fits'] = data[previousrun]['fits']

    # parse data into format for writing
    data_to_write = {} # only center position, for easy reading
    data_to_write_ext = {} # more information, for later plotting
    for run, val in data.items():
        run = int(run)
        data_to_write[run] = {}
        data_to_write_ext[run] = data[run]
        for varname, values in val['fits'].items():
            vartag = varname.split('_')[-1]
            data_to_write[run][vartag] = float(values[0])
    
    # write output data
    outputfile = os.path.join(args.outputdir, 'beamspot.json')
    with open(outputfile, 'w') as f:
        json.dump(data_to_write, f, indent=2)
    outputfile = os.path.join(args.outputdir, 'fitresults.json')
    with open(outputfile, 'w') as f:
        json.dump(data_to_write_ext, f, indent=2)
