import os
import sys
import glob
import time
import argparse
import subprocess

# local imports
thisdir = os.path.abspath(os.path.dirname(__file__))
topdir = os.path.abspath(os.path.join(thisdir, '..'))
sys.path.append(topdir)
from producetrees import read_samplelist


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, nargs='+',
      help='Input .root files, OR path to a .txt file listing input .root files (one per line)')
    parser.add_argument('-o', '--outputdir', required=True,
      help='Output directory')
    parser.add_argument('-n', '--nevents', type=int, default=-1,
      help='Number of events to process per file (default: all available events)')
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
    if not os.path.exists(args.outputdir): os.makedirs(args.outputdir)

    # loop over input files
    for input_file in input_files:

        # set output file name
        output_file = os.path.basename(input_file)
        output_file = os.path.join(args.outputdir, output_file)

        # make command
        cmd_stage1 = 'fccanalysis run analysis.py'
        if args.nevents > 0: cmd_stage1 += ' --nevents {}'.format(args.nevents)
        cmd_stage1 += ' --output {}'.format(output_file)
        cmd_stage1 += ' --files-list {}'.format(input_file)

        # create files storing stdout and stderr
        stdoutpath = output_file.replace('.root', '_stdout.txt')
        stderrpath = output_file.replace('.root', '_stderr.txt')
        stdout = open(stdoutpath, 'w')
        stderr = open(stderrpath, 'w')

        # run command
        print(f'Now running analysis...')
        print(cmd_stage1)
        start_time = time.time()
        subprocess.check_call(cmd_stage1, shell=True, stdout=stdout, stderr=stderr)
        end_time = time.time()
        stdout.write("Runtime: {:.3f}s.\n".format(end_time - start_time))
        stdout.close()
        stderr.close()
