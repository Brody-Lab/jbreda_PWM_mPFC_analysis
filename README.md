# jbreda_PWM_ephys_analysis

repository for analysis of PWM ephys data (as of 2020_1_6 this is for wireless mPFC tetrode recordings)

# TODO
* get preprocessed 0ing info
* `prep_protocol_info` currently a hard coded script, updated to flexible fx

## Assumptions

1. This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in kilosort. This is the best place to look if you want to understand how my files are structured. Additionally, I have taken & adjusted function from Tyler upload to bdata and Emily's [PWM](https://github.com/Brody-Lab/emilyjanedennis_PWManalysis/blob/master/find_wireless_sess.m) and [cscope](https://github.com/Brody-Lab/cscope) analyses repositories.

2. While you could implement some parts of this code pre or during sorting, I am implementing it after all my sorting has been completed for one rat. So, all the sessions I am interested in (from raw, to preprocessed, to sorted) a conveniently represented together in their respected directories. Thereby, most of these functions are intended to run *all* sessions for *one* rat as opposed to *some* sessions across *many* rats. Adjusting this should be feasible if you run a rat at a time.


## Session alignment & Kilosort/bdata extraction (matlab)
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

### utils (for matlab fxs)
(in no specific order)
1. `bundle_rename.m` is a script I used to rename all my sorted sessions from `session_name_preprocessinfo_firstbundle_forkilosort` to `session_name_bundle1_forkilosrt` to better align with Tyler's code

2. `mat_to_python_sess_info.ipynb` is a jupyter notebook that has code for important matlab structures and converting them to python dictionaries

3.  `ttl_match.mat` is an output from `find_wireless_sess_J` that I used in the `mat_to_python` notebook above and sometimes use it to check and remind myself what the output is :D

4. `readmda.m` is from mountainlab [git](https://github.com/flatironinstitute/mountainlab/tree/master/matlab/mdaio) & used to load .mda files into matlab

5. `kilosort_preprocess_mask_forcluster` is a function I edited to grab the masking data for a session (ie at what indices was the bundle zeroed out). As of 02/07/2021, I only ran on a single file and it will need to be adjusted to run on a batch of sessions.

### phy helpers

These functions are copied from [here](https://github.com/cortex-lab/spikes/tree/master/preprocessing/phyHelpers) & are primarily used to interact with kilosort output & load in matlab. Of note, I wrote an additional function `readClusterSpikeQualityCSV.m` because of how I sort my sessions. In short, **everything** that is a cell is marked in `cluster_group.tsv` as `good`, so I made in additional label called spike quaility (`sq`) saved as `cluster_sq.tsv`that informs whether a cell is multi or single. I've documented this well in `get_kspy_results_J.m` & even if you don't sort like this, you can use the function.


## Analysis (python)
```conda create -n PWM_ephys python=3.7 pip numpy matplotlib scipy scikit-learn h5py pyqt cython pillow

conda activate PWM_ephys
pip install spykes
pip install black
pip install seaborn
pip install jupyter

--- have but might not need ---
pip install spikeinterface (usually uses python=3.6)
pip install neo
```

### data_sdc_20190902_145404_exploratory_analysis.ipynb

This is the notebook I have started to write code in for plotting neural data aligned to behavior for a single session (20190902). I import the behavior & spike data extracted in the matlab functions above, wrangle them into a dataframe or dictionary, respectively, and then use a modified version of [Spykes](https://github.com/KordingLab/spykes) package to plot. It is a work in progress with current TODOs at the top. Most 'mature' functions can be found in `utils.py`. **The goal** of this is to eventually write a script that given the name of a session, can import behavior & spike info, and for a given set of events can get raster and psth info, plot it, and save out.

### utils.py

This contains functions I am current using in the notebook above & will eventually be used in the larger script

1.  `load_and_wrangle` which takes a path to behS and spkS and loads/wrangles them. Overwrite option if behavior has already been loaded & saved out as .csv Calls:
  a. `load_nested_mat` for dealing w/ the behS
  b. `load_behavior` loads behS into a dict
  c. `load_spks` loads spkS into a dict
  d. `make_beh_df` wrangles behavior dictionary into a df

2. `initiate_neurons` takes spike dictionary and creates `NeuroVis` object for each neuron

### spykes

This is a package form the Kording lab for making nice rasters/psths in python. I am currently using spyes.plot.NeuroVis. I have made some changes to their plotting functions to allow for axis to be assigned/titles to be made. Additionally, some of this code is outdated & should be re-written w/ seaborn integration eventually. To use their analysis/plotting tools, you need to pass spike times (in behavior time) into the `NeuroVis` class to create a usable object. You do not need to (and should not) event center your times to use this code; it will do that for you.
