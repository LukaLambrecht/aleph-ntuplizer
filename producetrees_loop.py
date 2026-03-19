# Looper over producetrees.py


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
    parser.add_argument('-o', '--outputdir', required=True,
      help='Output directory')
    parser.add_argument('--run-ntuplizer', default=False, action='store_true',
      help='Run per-jet ntuplizer stage (default: run only per-event stage)')
    parser.add_argument('--redirect-stdout', default=False, action='store_true',
      help='Redirect stdout and stderr from the terminal to txt files.')
    parser.add_argument('--do-clean', default=False, action='store_true',
      help='Remove intermediate output after running the ntuplizing stage.')
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
    if args.run_ntuplizer:
        cmd_compile = "g++ -o makentuples makentuples.cpp `root-config --cflags --libs` -Wall"
        print('Compiling makentuples...')
        subprocess.check_call(cmd_compile, shell = True, stdout=None, stderr=None)

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

        # find out if input file is simulation or data
        # (only used for output file naming)
        if 'QQB' in input_file: dtype = 'qqb'
        elif 'DATA' in input_file: dtype = 'data'
        else:
            msg = f'Data type of input file {input_file} not recognized.'
            raise Exception(msg)

        # make output file name
        outputfile = os.path.join(args.outputdir, f'output_{dtype}_{idx}.root')

        # make command to run
        cmd = 'python producetrees.py'
        cmd += f' -i {input_file}'
        cmd += f' -o {outputfile}'
        if args.run_ntuplizer:
            cmd += ' --run-ntuplizer'
            cmd += ' --no-compile'
        if args.redirect_stdout: cmd += ' --redirect_stdout'
        if args.do_clean: cmd += ' --do-clean'
        cmds.append(cmd)

    # run or submit commands
    if args.runmode == 'local':
        for cmd in cmds:
            print(cmd)
            os.system(cmd)
    elif args.runmode=='condor':
        env_script = os.path.abspath('../setup.sh')
        env_cmd = f'source {env_script}'
        ct.submitCommandsAsCondorCluster('cjob_producetrees', cmds,
          jobflavour='nextweek', conda_activate=env_cmd) 
