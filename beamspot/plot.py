import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

# global pyplot settings
plt.rc("text", usetex=True)
plt.rc("font", family="serif")


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', required=True,
      help='Input file, usually the output of the previous step (fit).')
    parser.add_argument('-o', '--outputdir', required=True,
      help='Output directory.')
    parser.add_argument('--sim', default=False, action='store_true',
      help='Run in simulation mode (only affects labels)')
    args = parser.parse_args()

    # make output directory
    if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

    # read input file
    with open(args.inputfile, 'r') as f:
        data = json.load(f)
    firstrun = list(data.keys())[0]
    varnames = list(data[firstrun]['fits'].keys())

    # make a summary figure
    fig, axs = plt.subplots(nrows=6, figsize=(12,12))
    colors = {'x': 'deepskyblue', 'y': 'royalblue', 'z': 'darkorchid'}
    for idx, varname in enumerate(varnames):
        ax1 = axs[2*idx]
        ax2 = axs[2*idx+1]
        means = np.array([data[run]['fits'][varname][0] for run in data.keys()])
        stds = np.array([data[run]['fits'][varname][1] for run in data.keys()])
        xax = np.arange(len(means))

        coord = varname.split('_')[-1]
        color = colors.get(coord, 'blue')

        # plot means and std
        ax1.fill_between(xax, means+stds, y2=means-stds, color=color, alpha=0.3)
        ax1.plot(xax, means, color=color, linewidth=2)
        mean = np.mean(means)
        if args.sim: ax1.axhline(0, linestyle='dashed', color='gray')
        else: ax1.axhline(mean, linestyle='dashed', color='gray')
        
        # plot std only
        ax2.fill_between(xax, stds, color=color, alpha=0.3)

        # plot aesthetics
        varlabel = f'Beamspot {coord}-coordinate'
        ax1.set_xticklabels([])
        ax1.set_xticks([])
        text = ax1.text(0.98, 0.95, varlabel + ' fitted center + width', ha='right', va='top',
                        transform=ax1.transAxes, fontsize=15)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))
        ax1.tick_params(labelsize=15)
        ax2.set_xticklabels([])
        ax2.set_xticks([])
        text = ax2.text(0.98, 0.95, varlabel + ' width', ha='right', va='top',
                        transform=ax2.transAxes, fontsize=15)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))
        ax2.grid(axis='y', which='both', linestyle='dashed', color='grey')
        ax2.set_ylim((0, ax2.get_ylim()[1]*1.2))
        ax2.tick_params(labelsize=15)

    # more plot aesthetics
    axs[5].set_xlabel('Run', fontsize=15)
    fig.subplots_adjust(wspace=0, hspace=0)

    # save figure
    fig.tight_layout()
    outputfile = os.path.join(args.outputdir, f'summary_beamspotfit.png')
    if args.sim: outputfile = outputfile.replace('.png', '_sim.png')
    fig.savefig(outputfile)
    fig.savefig(outputfile.replace('.png', '.pdf'))

    # make another summary figure
    fig, axs = plt.subplots(nrows=6, figsize=(12,12))
    colors = {'x': 'deepskyblue', 'y': 'royalblue', 'z': 'darkorchid'}
    for idx, varname in enumerate(varnames):
        ax1 = axs[2*idx]
        ax2 = axs[2*idx + 1]
        means = np.array([data[run]['fits'][varname][0] for run in data.keys()])
        uncs = np.array([data[run]['fits'][varname][2] for run in data.keys()])
        xax = np.arange(len(means))
        coord = varname.split('_')[-1]
        color = colors.get(coord, 'blue')

        # plot means and uncertainty
        ax1.fill_between(xax, means+uncs, y2=means-uncs, color=color, alpha=0.3)
        ax1.plot(xax, means, color=color)
        mean = np.mean(means)
        if args.sim: ax1.axhline(0, linestyle='dashed', color='gray')
        else: ax1.axhline(mean, linestyle='dashed', color='gray')

        # plot uncertainty only
        ax2.fill_between(xax, uncs, color=color, alpha=0.3)

        # plot aesthetics
        varlabel = 'Beamspot ' + r'$\bf{' + coord + '}$' + '-coordinate'
        ax1.set_xticklabels([])
        ax1.set_xticks([])
        text = ax1.text(0.98, 0.95, varlabel + ' fitted center',
          ha='right', va='top', transform=ax1.transAxes, fontsize=15)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))
        ax1.set_ylabel(f'{coord} [cm]', fontsize=15)
        ax1.tick_params(labelsize=15)
        
        ax2.set_xticklabels([])
        ax2.set_xticks([])
        text = ax2.text(0.98, 0.95, varlabel + ' uncertainty on fitted center',
          ha='right', va='top', transform=ax2.transAxes, fontsize=15)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))
        ax2.grid(axis='y', which='both', linestyle='dashed', color='grey')
        ax2.set_ylim((0, ax2.get_ylim()[1]*1.2))
        if coord=='x' or coord=='y': ax2.set_ylim((0, 0.002))
        elif coord=='z': ax2.set_ylim((0, 0.15))
        ax2.set_ylabel(f'{coord}-uncertainty [cm]', fontsize=15)
        ax2.tick_params(labelsize=15)

    # more plot aesthetics
    axs[5].set_xlabel('Run', fontsize=15)
    fig.subplots_adjust(wspace=0, hspace=0)

    # save figure
    fig.tight_layout()
    outputfile = os.path.join(args.outputdir, f'summary_beamspotcenter.png')
    if args.sim: outputfile = outputfile.replace('.png', '_sim.png')
    fig.savefig(outputfile)
    fig.savefig(outputfile.replace('.png', '.pdf'))
