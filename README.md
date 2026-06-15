# Ntuplizer for ALEPH data in EDM4HEP format


### Introduction
These scripts are intended to produce training and testing ntuples for a ParticleNet-like jet tagger study from EDM4HEP files (in particular the recently parsed ALEPH data).

This happens in two stages:
- `analysis.py`: Read the input data and calculate all properties of interest, stored in a per-event format.
- `makentuples.py`: Store all variables in per-jet ntuples ready for training a jet flavour classifier.


### Status
This repo is an attempt to make it run with the standard / default FCCAnalyses framework,
and undo the many custom changes that were applied to FCCAnalyses in an [earlier version](https://github.com/LukaLambrecht/FCCAnalyses/tree/aleph_ntuples) of this ntuplizer.

It is checked to be fully in sync with the latest ntuples produced with the earlier version (on 19 March 2026).
There are 2 branches `JetsConstituents_C` and `JetsConstituents_ct` with significant differences,
but that is ok becasue they were inconsistent in the earlier version and not used anywhere anyway.

All required modifications have been copied over to local analyzers,
such that the ntuplizer can run with the standard / default / vanilla FCCAnalyses framework.
Future work should go into trying to reduce the duplicated code and some hacks in the local analyzers,
potentially introducing small changes with respect to the earlier version
as long as they do not degrade the performance.
Another possibility is to backpropagate some of these changes and additions to the official FCCAnalyses framework.


### How to run
This code runs in a standard FCCAnalyses environment.
Activate it with `source /cvmfs/sw.hsf.org/key4hep/setup.sh` (e.g. on `lxplus`).
See more info [here](https://github.com/HEP-FCC/FCCAnalyses).

Then, run `python producetrees.py` with the appropriate options, for example:
```
python producetrees.py -i samplelists/samples_test.txt -o output_test/output.root -n 100
```

Note: these instructions are preliminary and the code might change extensively in the future.
You can always run `python producetrees.py -h` to see all available options.


### How to run on all available data

For simulation, use (for example, modify as needed):
```
python producetrees_loop.py -i samplelists/samples_alephsim.txt -o /eos/user/l/llambrec/aleph-data/ntuples-withdedx/eventlevel/mc -r condor --run-ntuplizer
```
This will run both stages in condor jobs (1 job per file, i.e. 70 jobs in total) and put the output files in the specified directory.
The output files of stage 1 (i.e. per-event storage) are called `output_qqb_<number>.root`.
The output files of stage 2 (i.e. flat ntuples with per-jet storage) are called `output_qqb_<number>_train.root`
and `output_qqb_<number>_test.root`.

Notes:
- As you can guess from the file names above,
an automatic split in training and testing set is performed, where 90\% of each input file is put in one file
(intended for training) and the remaining 10\% in another file (intended for testing).
This is not the most flexible solution and this behaviour might be changed at some point in the future (but it works fine for now).
- After running this, I usually go manually into the `ntuples` output folder, create `jetlevel/mc`,
and move the `*_train.root` and `*_test.root` from `eventlevel/mc` to `jetlevel/mc`.
This can of course be automated and/or made more flexible, just didn't get around to doing it yet.

For data, one can similarly use (for example):
```
python producetrees_loop.py -i samplelists/samples_alephdata.txt -o /eos/user/l/llambrec/aleph-data/ntuples-withdedx/eventlevel/data -r condor
```
which will similarly submit 1 job per file (i.e. 189 jobs in total).

Notes:
- While it is perfectly possible to add `--run-ntuplizer`, I usually don't,
as my subsequent analysis code is based on the per-event ntuples rather than the per-jet ntuples
(the latter ones are only used for training and testing the classifier).


### To do
- Make it more modular, allowing to disable some parts of the reconstruction with simple switches (for quicker testing). Currently it has to be done by defining and setting some global flag variables in `analysis.py`.
