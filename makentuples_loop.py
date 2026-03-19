# Looper over makentuples

# Use case: exceptional cases where I forgot to run the ntuplizing stage,
#           and just resubmitting the whole thing from scratch would take to long.


# external imports
import os
import sys
import six
import subprocess
import argparse

# local imports
thisdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(thisdir)
import tools.condortools as ct
from producetrees import read_samplelist


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, nargs='+',
      help='Input .root files, OR path to a .txt file listing input .root files (one per line)')
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

    # compile makentuples
    cmd_compile = "g++ -o makentuples makentuples.cpp `root-config --cflags --libs` -Wall"
    print('Compiling makentuples...')
    subprocess.check_call(cmd_compile, shell = True, stdout=None, stderr=None)

    # loop over input files
    cmds = []
    for idx, input_file in enumerate(input_files):

        # make ntuplizer commands 
        cmd_stagentuple_train = f'./makentuples {input_file} {input_file.replace(".root", "_train.root")}'
        cmd_stagentuple_train += f' 0 0.9'
        cmd_stagentuple_test = f'./makentuples {input_file} {input_file.replace(".root", "_test.root")}'
        cmd_stagentuple_test += f' 0.9 1'
        cmds.append(cmd_stagentuple_train)
        cmds.append(cmd_stagentuple_test)

    # run or submit commands
    if args.runmode == 'local':
        for cmd in cmds:
            print(cmd)
            os.system(cmd)
    elif args.runmode=='condor':
        env_script = os.path.abspath('../setup.sh')
        env_cmd = f'source {env_script}'
        ct.submitCommandsAsCondorCluster('cjob_producetrees', cmds,
          jobflavour='workday', conda_activate=env_cmd) 
