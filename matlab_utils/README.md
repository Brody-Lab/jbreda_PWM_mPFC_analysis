1. `bundle_rename.m` is a script I used to rename all my sorted sessions from `session_name_preprocessinfo_firstbundle_forkilosort` to `session_name_bundle1_forkilosrt` to better align with Tyler's code

2.  `ttl_match.mat` is an output from `find_wireless_sess_J` that I used in the `mat_to_python` notebook above and sometimes use it to check and remind myself what the output is :D

3. `readmda.m` is from mountainlab [git](https://github.com/flatironinstitute/mountainlab/tree/master/matlab/mdaio) & used to load .mda files into matlab

4. `kilosort_preprocess_mask_forcluster` is a function I edited to grab the masking data for a session (ie at what indices was the bundle zeroed out). As of 06/01/2020, this is not used.

5. `phy helpers` These functions are copied from [here](https://github.com/cortex-lab/spikes/tree/master/preprocessing/phyHelpers) & are primarily used to interact with kilosort output & load in matlab. Of note, I wrote an additional function `readClusterSpikeQualityCSV.m` because of how I sort my sessions. In short, **everything** that is a cell is marked in `cluster_group.tsv` as `good`, so I made in additional label called spike quaility (`sq`) saved as `cluster_sq.tsv`that informs whether a cell is multi or single. I've documented this well in `get_kspy_results_J.m` & even if you don't sort like this, you can use the function.
