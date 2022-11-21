## Session alignment & Kilosort/bdata extraction (matlab)

Functions to be run in order
----

**`find_wireless_sess_J.m` :**
---
takes a session, finds the corresponding behavioral session in the`bdata` SQL table and outputs information to use for alignment

**specific details:**

This function compares TTL ITI from the headstage (trodes) and computer (finite state machine) and finds the behavioral session within a 3 day window of the trodes date with the smallest residual via linear regression. Trodes ttl `.dio` file is read in via `readTrodesExtractedDataFile.m` 

This function specific to the PWM task and uses the wave `TrigScope` for alignment. Using the output from this regression of headstage ttl vs fsm ttls, this function saves a function to convert from spike times to behavior times in `ttl_match.mat`. 

Additional helpful things saved are: session time, sessionid in bdata, useful paths for alignment, and a figure of ttls from both systems for the chosen session under `ttl_alignment.jpeg`, among others.

If you want to run multiple sessions through this function --> `find_all_wireless_sess_J.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `find_wireless_sess_J`


**`get_ksphy_results_J.m` :** 
---
takes a session, loads the `kilosort` information for each bundle in that session, loads the `.mda` file for tetrodes with cells, and grabs & saves spike information for each cell in the session.

**specific details**:

This function uses uses `phyHelpers` and  from output from `fine_wirless_sess_J` to find & load relevant `kilosort` information for the session (e.g. cluster id, single/multi unit). For each cell, the strongest tetrode .mda file is loaded & filtered. Then, spike times, spike indices, average waveforms are saved out along with an average waveform. 

The functions saves into two different structs with very similar information, but different shapes. `spkS` saves as `ksphy_clusters_forbdata.mat` is specific to Tyler's code & uploading to bdata. `PWMspkS` saves as `ksphy_clusters_foranalysis.mat` is a struct I wrote for analysis in python.

If you want to run multiple sessiosn through this function --> `find_all_ksphy_results.m` takes a directory with session names (ie the directory that has all my sorted sessions) and runs them through `get_ksphy_results`

**`prep_protocol_info.m`:**
---
takes a session, loads the behavior session performed w/ ephys & extracts protocol information to be used for spike time analysis.

**specific details**:

This function takes a behavior session, runs `get_ksphy_results_J` to get behavior alignment info (could also run `find_wireless_sess_J` instead if you want), and loads the correct behavior session from `bdata`. 

It then pulls out a variety of relevant trial information for PWM (eg n_done_trials, sound_delay_time), and saves them into a struct `behS` as `protocol_info.mat`. Along with the functions above, this was done for the sake of export into python & load into dictionaries.

### Matlab_utils

These are utils used for the steps above that come from different code libraries/labs. They have there own README with detailed information, check it out **here**