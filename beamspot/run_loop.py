# Looper over run.py


# external imports
import os
import sys
import six
import argparse

# local imports
thisdir = os.path.abspath(os.path.dirname(__file__))
topdir = os.path.abspath(os.path.join(thisdir, '..'))
sys.path.append(topdir)
import tools.condortools as ct
from producetrees import read_samplelist


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, nargs='+',
      help='Input .root files, OR path to a .txt file listing input .root files (one per line)')
    parser.add_argument('-o', '--outputdir', required=True,
      help='Output directory')
    parser.add_argument('-r', '--runmode', default='local', choices=['local', 'condor'])
    args = parser.parse_args()

    # find input files
    input_files = []
    for el in args.input:
        if el.endswith('.root'): input_files.append(el)
        elif el.endswith('.txt'):
            input_files += read_samplelist(el)
    print(f'Found following input files ({len(input_files)}):')
    for f in input_files: print(f'  - {f}')

    # make output directory
    if os.path.exists(args.outputdir):
        msg = f'Output directory {args.outputdir} already exists; clean? (y/n)'
        print(msg)
        go = six.moves.input()
        if go!='y': sys.exit()
        os.system(f'rm {os.path.join(args.outputdir, "*")}')
    else: os.makedirs(args.outputdir)

    # loop over input files
    cmds = []
    for idx, input_file in enumerate(input_files):

        # make command to run
        cmd = 'python run.py'
        cmd += f' -i {input_file}'
        cmd += f' -o {args.outputdir}'
        cmds.append(cmd)

    # run or submit commands
    if args.runmode == 'local':
        for cmd in cmds:
            print(cmd)
            os.system(cmd)
    elif args.runmode=='condor':
        env_script = os.path.abspath('../../setup.sh')
        env_cmd = f'source {env_script}'
        ct.submitCommandsAsCondorCluster('cjob_producetrees', cmds,
          jobflavour='workday', conda_activate=env_cmd) 
