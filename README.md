# jbreda_PWM_mPFC_analysis

repository for analysis of rat mPFC PWM ephys data for tetrode recordings

## Assumptions

1. This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in kilosort. This is the best place to look if you want to understand how my files are structured. Additionally, I have taken & adjusted function from Tyler's upload to bdata and Emily's [PWM](https://github.com/Brody-Lab/emilyjanedennis_PWManalysis/blob/master/find_wireless_sess.m) and [cscope](https://github.com/Brody-Lab/cscope) analyses repositories.

2. While you could implement some parts of this code pre or during sorting, I am implementing it after all my sorting has been completed for one rat. So, all the sessions I am interested in (from raw, to preprocessed, to sorted) a conveniently represented together in their respected directories. Thereby, most of these functions are intended to run *all* sessions for *one* rat as opposed to *some* sessions across *many* rats. Adjusting this should be feasible if you run a rat at a time.


## Session alignment & Kilosort/bdata extraction (matlab)
(in order)
1. `find_wireless_sess_J.m` takes a session, finds the corresponding behavioral session in bdata and output information to use for alignment

  **specific details**:

  This function compares TTL ITI from the headstage (trodes) and computer (finite state machine) and finds the behavioral session within a 3 day window of the trodes date with the smallest residual via linear regression. Trodes ttl `.dio` file is read in via `readTrodesExtractedDataFile.m` This function specific to the PWM task and uses the wave `TrigScope` for alignment. Using the output from this regression of headstage ttl vs fsm ttls, this function saves a function to convert from spike times to behavior times in `ttl_match.mat`. Additional helpful things saved are: session time, sessionid in bdata, useful paths for alignment, and a figure of ttls from both systems for the chosen session under `ttl_alignment.jpeg`, among others.

  If you want to run multiple sessions through this function --> `find_all_wireless_sess_J.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `find_wireless_sess_J`


2. `get_ksphy_results_J.m` takes a session, loads the kilosort information for each bundle in that session, loads the .mda file for tetrodes with cells, and grabs & saves spike information for each cell in the session.

  **specific details**:

  This function uses uses `phyHelpers` and  from output from `fine_wirless_sess_J` to find & load relevant kilosort information for the session (e.g. cluster id, single/multi unit). For each cell, the strongest tetrode .mda file is loaded & filtered. Then, spike times, spike indices, average waveforms are saved out along with an average waveform. The functions saves into two different structs with very similar information, but different shapes. `spkS` saves as `ksphy_clusters_forbdata.mat` is specific to Tyler's code & uploading to bdata. `PWMspkS` saves as `ksphy_clusters_foranalysis.mat` is a struct I wrote for analysis in python.

  If you want to run multiple sessiosn through this function --> `find_all_ksphy_results.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `get_ksphy_results`

3. `prep_protocol_info.m` takes a session, loads the behavior session performed w/ ephys & extracts protocol information to be used for spike time analysis.

  **specific details**:

  This function takes a behavior session, runs `get_ksphy_results_J` to get behavior alignment info (could also run `find_wiresless_sess_J` instead if you want), and loads the correct behavior session from bdata. Then, it pulls out a variety of relevant trial information for PWM (eg n_done_trials, sound_delay_time), and saves them into a struct `behS` as `protocol_info.mat`. Along with the functions above, this was done for the sake of export into python & load into dictionaries.

### Matlab_utils

These are utils I use for the steps above that come from different code libraries/labs. They have there own README with detailed information, look there if you'd like.

## Analysis (python)

Environment information. Note the jedi downgrade is for a tab-complete bug

```
- conda create -n PWM_ephys python=3.7 pip numpy matplotlib scipy scikit-learn h5py pyqt cython pillow black seaborn jupyter statsmodels pydove
conda activate PWM_ephys
pip3 install jedi==0.17.2

```

### Current Status

This repo is going to sit for a few months while we collect data from more animals and then will either be forked into a new repo or revived.

### Tutorial

To see how I was using the code as of June 2021, see the `2021_Tutorial.ipynb`. Functions are either located in `io_utils.py` for importing/wrangling or `plotting_utils.py` for figuring making/analysis.

This code primarily focuses on analyzing the delay period. It initially was built off of the [Spykes](https://github.com/jess-breda/spykes) library, but is now all written from scratch. There were a lot of issues with noise and kilosort fully masking some of the data. With smoothing, this has become less of an issue, and all the data in the tutorial is unmasked.

### Future Upgrades

Making the aligned spiking times for each neuron/event into a python object that can easily be filtered for whatever conditions your interested in (e.g. hit only, go left) would be great. Right now, the code needs to be run for each condition you want (see tutorial for example).

### Final Notes  

The `scraps` folder contains notebooks/other folders that I used at one point and I'm 99% sure I will not need again, but I am waiting to throw them out until I start a new round of analyses.

The `analysis_for_presentations` folder has notebooks I used to generate figures for lab meeting & Carlos meetings
