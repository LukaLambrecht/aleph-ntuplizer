import os
import sys
import glob
import time
import argparse
import subprocess
import multiprocessing as mp
from subprocess import Popen,PIPE
from datetime import date


def read_samplelist(samplelist):
    '''
    Read a sample list.
    Input arguments:
    - samplelist: path to txt file listing samples.
      every line in the sample list is assumed to be a path to a root file
      (potentially containing unix-style wildcards).
    '''
    # read each line from the sample list
    with open(samplelist, 'r') as f:
        lines = f.readlines()
    # cleaning and filtering
    lines = [l.strip(' \t\n') for l in lines]
    lines = [l for l in lines if (l.endswith('.root') and not l.startswith('#'))]
    # expand wildcards
    input_files = sum([glob.glob(l) for l in lines], [])
    return input_files


def run_ntuplizer(cmd_stagentuple_train, cmd_stagentuple_test, f_stdout, f_stderr):
    '''Internal helper function to run the ntuplizer commands'''
    start_time = time.time()
    subprocess.check_call(cmd_stagentuple_train, shell=True, stdout=f_stdout, stderr=f_stderr)
    subprocess.check_call(cmd_stagentuple_test, shell=True, stdout=f_stdout, stderr=f_stderr)
    end_time = time.time()
    msg = 'Ntuplizing stage runtime: {:.3f}s.'.format(end_time - end_time)
    if f_stdout is None: print(msg)
    else: f_stdout.write(msg+'\n')


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, nargs='+',
      help='Input .root files, OR path to a .txt file listing input .root files (one per line)')
    parser.add_argument('-o', '--outputfile', required=True,
      help='Output .root file')
    parser.add_argument('-n', '--nevents', type=int, default=-1,
      help='Number of events to process (default: all available events)')
      # note: when multiple files are provided, this number seems to be the total number of events,
      #       not the number of events per file!
      #       so some files may not be read at all if nevents is reached in an earlier file.
    parser.add_argument('--run-ntuplizer', default=False, action='store_true',
      help='Run per-jet ntuplizer stage (default: run only per-event stage)')
    parser.add_argument('--training-frac', type=float, default=0.9,
      help='Fraction of events to put in training ntuple (default: 0.9)'
          +' (rest will go in testing ntuple)')
      # todo: find out if shuffling would be needed here.
    parser.add_argument('--redirect-stdout', default=False, action='store_true',
      help='Redirect stdout and stderr from the terminal to txt files.')
    parser.add_argument('--do-clean', default=False, action='store_true',
      help='Remove intermediate output after running the ntuplizing stage.')
    args = parser.parse_args()

    # find input files
    input_files = []
    for el in args.input:
        if el.endswith('.root'): input_files.append(el)
        elif el.endswith('.txt'):
            input_files += read_samplelist(el)
    print(f'Found following input files ({len(input_files)}):')
    for f in input_files: print(f'  - {f}')

    # check output file
    # (naming depends on output file ending in '.root', maybe make more robust/flexible later)
    if not args.outputfile.endswith('.root'):
        msg = 'Output file must end with ".root"'
        raise Exception(msg)

    # make output directory
    outputdir = os.path.dirname(args.outputfile)
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    # make command for stage 1
    tempfile = args.outputfile
    cmd_stage1 = 'fccanalysis run analysis.py'
    if args.nevents > 0: cmd_stage1 += ' --nevents {}'.format(args.nevents)
    cmd_stage1 += ' --output {}'.format(tempfile)
    cmd_stage1 += ' --files-list {}'.format(' '.join(input_files))

    # make command for stage 2 
    cmd_stagentuple = f'python makentuples.py -i {tempfile}'

    # create files storing stdout and stderr
    stdoutpath = args.outputfile.replace('.root', '_stdout.txt')
    stderrpath = args.outputfile.replace('.root', '_stderr.txt')
    stdout = None
    stderr = None
    if args.redirect_stdout:
        stdout = open(stdoutpath, 'w')
        stderr = open(stderrpath, 'w')

    # run stage 1
    msg = f'Now running stage 1:\n{cmd_stage1}'
    if stdout is None: print(msg)
    else: stdout.write(msg+'\n')
    start_time = time.time()
    subprocess.check_call(cmd_stage1, shell=True, stdout=stdout, stderr=stderr)
    end_time = time.time()
    msg = 'Stage 1 runtime: {:.3f}s.'.format(end_time - start_time)
    if stdout is None: print(msg)
    else: stdout.write(msg+'\n')

    # run stage 2
    if args.run_ntuplizer:
        threads = []
        cmd_stagentuple_train = cmd_stagentuple
        cmd_stagentuple_train += ' -o ' + args.outputfile.replace('.root', '_train.root')
        cmd_stagentuple_train += f' --entry_start 0 --entry_stop {args.training_frac}'
        cmd_stagentuple_test = cmd_stagentuple
        cmd_stagentuple_test += ' -o ' + args.outputfile.replace('.root', '_test.root')
        cmd_stagentuple_test += f' --entry_start {args.training_frac} --entry_stop 1'
        msg = f'Now running stage 2:\n{cmd_stagentuple_train}\n{cmd_stagentuple_test}'
        if stdout is None: print(msg)
        else: stdout.write(msg+'\n')
        thread = mp.Process(target=run_ntuplizer, args=(cmd_stagentuple_train, cmd_stagentuple_test, stdout, stderr))
        thread.start()
        threads.append(thread)
        
        for proc in threads: proc.join()

        # do cleaning: remove stage 1 output
        if args.do_clean:
            print('Removing output from stage 1 for cleaning...')
            os.system(f'rm {tempfile}')

    if stdout is not None: stdout.close()
    if stderr is not None: stderr.close()

    # do cleaning: remove stdout and stderr log files
    if args.do_clean:
        if stdout is not None: os.system(f'rm {stdoutpath}')
        if stderr is not None: os.system(f'rm {stderrpath}')
