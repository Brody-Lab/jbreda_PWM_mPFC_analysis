# jbreda_PWM_ephys_analysis

repository for analysis of PWM ephys data (as of 2020_1_6 this is for wireless mPFC tetrode recordings)

# TODO
* get preprocessed 0ing info
* `prep_protocol_info` currently a hard coded script, updated to flexible fx

## Assumptions

1. This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in kilosort. This is the best place to look if you want to understand how my files are structured. Additionally, I have taken & adjusted function from Tyler upload to bdata and Emily's [PWM](https://github.com/Brody-Lab/emilyjanedennis_PWManalysis/blob/master/find_wireless_sess.m) and [cscope](https://github.com/Brody-Lab/cscope) analyses repositories.

2. While you could implement some parts of this code pre or during sorting, I am implementing it after all my sorting has been completed for one rat. So, all the sessions I am interested in (from raw, to preprocessed, to sorted) a conveniently represented together in their respected directories. Thereby, most of these functions are intended to run *all* sessions for *one* rat as opposed to *some* sessions across *many* rats. Note that adjusting this should be quite easy.


## Session alignment
(in order)
1. `find_wireless_sess_J.m` takes a session, finds the corresponding behavioral session in bdata and output information to use for alignment

  **specific details**:

  This function compares TTL ITI from the headstage (trodes) and computer (finite state machine) and finds the behavioral session within a 3 day window of the trodes date with the smallest residual via linear regression. Trodes ttl `.dio` file is read in via `readTrodesExtractedDataFile.m` This function specific to the PWM task and uses the wave `TrigScope` for alignment. Using the output from this regression of headstage ttl vs fsm ttls, this function saves a function to convert from spike times to behavior times in `ttl_match.mat`. Additional helpful things saved are: session time, sessionid in bdata, useful paths for alignment, and a figure of ttls from both systems for the chosen session under `ttl_alignment.jpeg`, among others.

  If you want to run multiple sessions through this function --> `find_all_wireless_sess.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `find_wireless_sess_J`


2. `get_ksphy_results_J.m` takes a session, loads the kilosort information for each bundle in that session, loads the .mda file for tetrodes with cells, and grabs & saves spike information for each cell in the session.

  **specific details**:

  This function uses uses `phyHelpers` and  from output from `fine_wirless_sess_J` to find & load relevant kilosort information for the session (e.g. cluster id, single/multi unit). For each cell, the strongest tetrode .mda file is loaded & filtered. Then, spike times, spike indices, average waveforms are saved out along with an average waveform. The functions saves into two different structs with very similar information, but different shapes. `spkS` saves as `ksphy_clusters_forbdata.mat` is specific to Tyler's code & uploading to bdata. `PWMspkS` saves as `ksphy_clusters_foranalysis.mat` is a struct I wrote for analysis in python.

  If you want to run multiple sessiosn through this function --> `find_all_ksphy_results.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `get_ksphy_results`

3. `prep_protocol_info.m` takes a session, loads the behavior session performed w/ ephys & extracts protocol information to be used for spike time analysis.

  **specific details**:

  This function takes a behavior session, runs `get_ksphy_results_J` to get behavior alignment info (could also run `find_wiresless_sess_J` instead if you want), and loads the correct behavior session from bdata. Then, it pulls out a variety of relevant trial information for PWM (eg n_done_trials, sound_delay_time), and saves them into a struct `behS` as `protocol_info.mat`. Along with the functions above, this was done for the sake of export into python & load into dictionaries.


## python
```conda create -n PWM_ephys python=3.7 pip numpy matplotlib scipy scikit-learn h5py pyqt cython pillow

conda activate PWM_ephys
pip install neo
pip install spykes
pip install spikeinterface (usually uses python=3.6)
pip install neo
pip install black
pip install seaborn
```

* working on exploratory analysis for 20190902 in notebook, will update when done

## utils
(in no specific order)
1. `bundle_rename.m` is a script I used to rename all my sorted sessions from `session_name_preprocessinfo_firstbundle_forkilosort` to `session_name_bundle1_forkilosrt` to better align with Tyler's code

2. `mat_to_python_sess_info.ipynb` is a jupyter notebook that has code for important matlab structures and converting them to python dictionaries

3.  `ttl_match.mat` is an output from `find_wireless_sess_J` that I used in the `mat_to_python` notebook above and sometimes use it to check and remind myself what the output is :D

4. `readmda.m` is from mountainlab [git](https://github.com/flatironinstitute/mountainlab/tree/master/matlab/mdaio) & used to load .mda files into matlab

5. `kilosort_preprocess_mask_forcluster`

## phy helpers

These functions are copied from [here](https://github.com/cortex-lab/spikes/tree/master/preprocessing/phyHelpers) & are primarily used to interact with kilosort output & load in matlab. Of note, I wrote an additional function `readClusterSpikeQualityCSV.m` because of how I sort my sessions. In short, **everything** that is a cell is marked in `cluster_group.tsv` as `good`, so I made in additional label called spike quaility (`sq`) saved as `cluster_sq.tsv`that informs whether a cell is multi or single. I've documented this well in `get_kspy_results_J.m` & even if you don't sort like this, you can use the function.
