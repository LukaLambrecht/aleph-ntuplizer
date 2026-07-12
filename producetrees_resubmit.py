# Resubmit failed jobs from a previous `producetrees_loop.py -r condor` run.
#
# Deliberately does not touch producetrees_loop.py or rely on the condor submit
# description file it writes (tools/condortools.py silently overwrites that file
# in full on every new submission that reuses the same job name, so it cannot be
# trusted as a record of an earlier submission). Instead, this reads the actual
# per-job stdout log files condor already writes (cjob_..._out_<ClusterId>_<ProcId>),
# which are inherently protected against collision by their ClusterId_ProcId suffix,
# and which already contain the exact input/output file(s) producetrees.py used
# (printed via its own "Now running stage 1: fccanalysis run analysis.py --output
# ... --files-list ..." line).
#
# This is a deliberately simple, standalone stopgap -- not wired into
# producetrees_loop.py at all. Can be revisited more cleanly later.
#
# Usage:
#   python producetrees_resubmit.py [--run-ntuplizer] [--do-clean] [-r local|condor]
#
# WARNING: only run this once the original submission has fully finished (e.g.
# confirmed via condor_q). A missing/invalid output can mean either a failed job
# or one still queued/running -- this script cannot tell the difference, and
# resubmitting a still-running job creates a duplicate job writing to the same
# output file.


import os
import re
import sys
import glob
import six
import uproot
import argparse

# local imports
thisdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(thisdir)
import tools.condortools as ct


def is_valid_root_file(filepath, treename='events', min_entries=1):
    '''
    Check whether filepath is a readable ROOT file containing a tree named
    treename with at least min_entries entries.
    Returns False for a missing file, a file uproot cannot open (e.g. truncated
    or corrupted, as can happen from a transient filesystem glitch), a file
    missing the expected tree, or too few entries in it.
    '''
    if not os.path.exists(filepath): return False
    try:
        with uproot.open(filepath) as f:
            if treename not in f: return False
            return f[treename].num_entries >= min_entries
    except Exception:
        return False


def parse_job_log(logfile):
    '''
    Parse a cjob_..._out_<ClusterId>_<ProcId> stdout log file (as produced by
    producetrees.py running under condor) for the input file(s) and output file
    it was given, from its own
    "fccanalysis run analysis.py [--nevents N] --output <output> --files-list <input(s)>"
    print statement (see producetrees.py), and additionally for the number of
    events fccanalysis itself reports having processed ("Total events processed:"
    in its own SUMMARY block).
    Returns (inputfiles, outputfile, nevents_processed); any of these is None if
    not found (e.g. the job never got far enough to print it).
    Note on nevents_processed: this is needed because a file can be *openable*
    with a non-empty tree and still be silently truncated -- ROOT's automatic
    file-recovery mechanism (visible as "TFile::Recover" messages in these same
    log files) can patch up a file that a filesystem glitch interrupted mid-write
    into something that opens fine but contains far fewer events than were
    actually processed, so a plain "does it open, is it non-empty" check is not
    enough on its own to catch that case.
    '''
    with open(logfile, 'r') as f:
        content = f.read()
    match = re.search(r'--output\s+(\S+)\s+--files-list\s+(.+)', content)
    if not match:
        return None, None, None
    outputfile = match.group(1)
    inputfiles = match.group(2).split()
    nevents_processed = None
    nmatch = re.search(r'Total events processed:\s*([\d,]+)', content)
    if nmatch:
        nevents_processed = int(nmatch.group(1).replace(',', ''))
    return inputfiles, outputfile, nevents_processed


def find_job_suffixes(logdir):
    '''
    Find all distinct <ClusterId>_<ProcId> suffixes among cjob_*_out_*/cjob_*_err_*
    files in logdir (looking at both in case a job is for some reason missing one
    of the two).
    '''
    suffixes = set()
    for tag in ('out', 'err'):
        for fname in glob.glob(os.path.join(logdir, f'cjob_*_{tag}_*')):
            base = os.path.basename(fname)
            idx = base.rfind(f'_{tag}_')
            if idx < 0: continue
            suffixes.add(base[idx+len(f'_{tag}_'):])
    return sorted(suffixes)


if __name__ == '__main__':

    # read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-ntuplizer', default=False, action='store_true',
      help='Run per-jet ntuplizer stage (default: run only per-event stage);'
          +' same meaning as in producetrees_loop.py. Determines which output'
          +' file(s) are checked/resubmitted with (train/test vs. plain per-event).')
    parser.add_argument('--do-clean', default=False, action='store_true',
      help='Remove intermediate output after running the ntuplizing stage;'
          +' same meaning as in producetrees_loop.py.')
    parser.add_argument('-r', '--runmode', default='local', choices=['local', 'condor'])
    parser.add_argument('--jobflavour', default='workday')
    parser.add_argument('--logdir', default='.',
      help='Directory to scan for cjob_* log files (default: current directory).')
    parser.add_argument('--jobname', default='cjob_producetrees_resubmit',
      help='Base name for the new condor job files, if -r condor.')
    args = parser.parse_args()

    # find and parse all job log files
    suffixes = find_job_suffixes(args.logdir)
    print(f'Found {len(suffixes)} job log file(s) (by ClusterId_ProcId) in {args.logdir}.')

    njobs = 0
    nvalid = 0
    cmds = []
    unparsed = []
    for suffix in suffixes:
        candidates = glob.glob(os.path.join(args.logdir, f'cjob_*_out_{suffix}'))
        if len(candidates) == 0:
            unparsed.append(suffix)
            continue
        inputfiles, outputfile, nevents_processed = parse_job_log(candidates[0])
        if outputfile is None:
            unparsed.append(suffix)
            continue
        njobs += 1

        # determine expected output(s), same logic as producetrees.py itself:
        # if --run-ntuplizer, the train/test files matter (tree 'tree'), not the
        # intermediate per-event file (which --do-clean may have removed by design);
        # otherwise the plain per-event file (tree 'events'). For the plain
        # per-event case, also cross-check against the number of events the log
        # says fccanalysis actually processed, not just "non-empty" -- see the
        # note in parse_job_log() on why a truncated-but-recovered file can
        # otherwise look valid.
        if args.run_ntuplizer:
            expected = [(outputfile.replace('.root', '_train.root'), 'tree', 1),
                        (outputfile.replace('.root', '_test.root'), 'tree', 1)]
        else:
            min_entries = nevents_processed if nevents_processed else 1
            expected = [(outputfile, 'events', min_entries)]

        if all(is_valid_root_file(path, tname, minent) for path, tname, minent in expected):
            nvalid += 1
            continue

        # build the resubmission command
        cmd = 'python producetrees.py'
        cmd += ' -i ' + ' '.join(inputfiles)
        cmd += f' -o {outputfile}'
        if args.run_ntuplizer: cmd += ' --run-ntuplizer'
        if args.do_clean: cmd += ' --do-clean'
        cmds.append(cmd)

    print(f'{njobs} job(s) successfully parsed; {len(unparsed)} could not be parsed'
         +' (not resubmitted automatically -- check these by hand):')
    for s in unparsed: print(f'  - {s}')
    print(f'{nvalid} already have valid output; {len(cmds)} will be (re)submitted:')
    for cmd in cmds: print(f'  - {cmd}')

    if len(cmds) == 0:
        print('Nothing to resubmit.')
        sys.exit()

    # important safety note: a job with no valid output yet could either have
    # failed, or simply not be finished (or not even started) yet -- this script
    # has no way to tell the difference (e.g. no use of condor_q here). Resubmitting
    # a job that is still queued or running would create a duplicate, competing
    # condor job writing to the same output path.
    print()
    print('WARNING: a missing/invalid output can mean either that the job failed,')
    print('OR that it is simply still queued/running and has not produced output yet.')
    print('Only proceed once you have confirmed (e.g. via condor_q) that the ORIGINAL')
    print('submission has fully finished -- resubmitting still-running jobs creates')
    print('duplicate jobs writing to the same output files.')
    print()
    print('Continue? (y/n)')
    go = six.moves.input()
    if not go == 'y': sys.exit()

    # run or submit commands
    if args.runmode == 'local':
        for cmd in cmds:
            print(cmd)
            os.system(cmd)
    elif args.runmode == 'condor':
        env_script = os.path.abspath('setup.sh')
        env_cmd = f'source {env_script}'
        ct.submitCommandsAsCondorCluster(args.jobname, cmds,
          jobflavour=args.jobflavour, conda_activate=env_cmd)
