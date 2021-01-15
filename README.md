# jbreda_PWM_ephys_analysis

repository for analysis of PWM ephys data (as of 2020_1_6 this is for wireless mPFC tetrode recordings)

# TODO
* get preprocessed 0ing info
* adjust get_ksphy_results to deal w/ sq
* assess if filtering in get_ksphy_results is appropriate for these sessions

## Assumptions

1. This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in kilosort. This is the best place to look if you want to understand how my files are structured. Additionally, I have taken & adjusted function from Tyler upload to bdata and Emily's [PWM](https://github.com/Brody-Lab/emilyjanedennis_PWManalysis/blob/master/find_wireless_sess.m) and [cscope](https://github.com/Brody-Lab/cscope) analyses repositories.

2. While you could implement some parts of this code pre or during sorting, I am implementing it after all my sorting has been completed for one rat. So, all the sessions I am interested in (from raw, to preprocessed, to sorted) a conveniently represented together in their respected directories. Thereby, most of these functions are intended to run *all* sessions for *one* rat as opposed to *some* sessions across *many* rats. Note that adjusting this should be quite easy.


## Session alignment
(in order)
1. `find_wireless_sess_J.m` takes a session, finds the corresponding behavioral session in bdata and output information to use for alignment

  **specific details**:

  This function compares TTL ITI from the headstage (trodes) and computer (finite state machine) and finds the behavioral session within a 3 day window of the trodes date with the smallest residual via linear regression. Trodes ttl `.dio` file is read in via `readTrodesExtractedDataFile.m` This function specific to the PWM task and uses the wave `TrigScope` for alignment. Using the output from this regression of headstage ttl vs fsm ttls, this function saves a function to convert from spike times to behavior times in `ttl_match.mat`. Additional helpful things saved are: session time, sessionid in bdata, useful paths for alignment, and a figure of ttls from both systems for the chosen session under `ttl_alignment.jpeg`, among others.

2. `find_all_wireless_sess.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `find_wireless_sess_J`


3. `get_ksphy_results.m`
***WORKING HERE***
- got it to run on a session on 1/15/2020
- using a lot of `phyHelpers` from [here](https://github.com/cortex-lab/spikes/tree/master/preprocessing/phyHelpers) and wrote my own `readClusterSpikeQualityCSV.m`
- issues with sq sometimes having a 0 output and that messes up single/multi assignments further down
- Line 160 on needs a closer look
  - is do I want to do that filter..?


## utils
(in no specific order)
1. `bundle_rename.m` is a script I used to rename all my sorted sessions from `session_name_preprocessinfo_firstbundle_forkilosort` to `session_name_bundle1_forkilosrt` to better align with Tyler's code

2. `mat_to_python_sess_info.ipynb` is a jupyter notebook that has code for important matlab structures and converting them to python dictionaries

3.  `ttl_match.mat` is an output from `find_wireless_sess_J` that I used in the `mat_to_python` notebook above and sometimes use it to check and remind myself what the output is :D
