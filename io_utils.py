#!/usr/bin/env python

"""
Utils to be used for importing & wrangling spiking & behavior data curated in
kilosort & taken from Brodylab bdatabase
"""

# libraries
import sys; sys.path.insert(0, '..') # if you don't find it here, look one above
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.io as spio

## ===Importing .mat structs & extracting ===
def load_nested_mat(filename):
    """
    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects. From https://stackoverflow.com/questions
    /48970785/complex-matlab-struct-mat-file-read-by-python
    """

    def _check_vars(d):
        """
        Checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in d:
            if isinstance(d[key], spio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
            elif isinstance(d[key], np.ndarray):
                d[key] = _toarray(d[key])
        return d

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _toarray(elem)
            else:
                d[strg] = elem
        return d

    def _toarray(ndarray):
        """
        A recursive function which constructs ndarray from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        if ndarray.dtype != 'float64':
            elem_list = []
            for sub_elem in ndarray:
                if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                    elem_list.append(_todict(sub_elem))
                elif isinstance(sub_elem, np.ndarray):
                    elem_list.append(_toarray(sub_elem))
                else:
                    elem_list.append(sub_elem)
            return np.array(elem_list)
        else:
            return ndarray

    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_vars(data)


def load_and_wrangle(beh_path, spks_path, overwrite):
    """
    This function loads behavior and spike data from a single session given path
    information. If behavior data has already been loaded & wrangled, will load
    the dataframe instead of creating a new one.

    inputs
    ------
    beh_path  : string, path to .mat file with behavior info
    spks_path : string, path to .mat file with spks info
    overwrite (optional, False) : wether or not to overwrite previous dataframe
    or load up if already made

    returns
    ------
    beh_df    :  df (ntrials x items), tidy data frame with behavior information
                & some relabeling
    spks_info : ndarray, with spk info from .mat file
    """

    # --Spikes--- (eventually can load/in out with pickle if needed)
    spks_info = load_spks(spks_path)

    # --Behavior--
    # check if the df has already been created, then either load or wrangle
    session_dir = os.path.dirname(beh_path)

    #this is broken & I need to fix
    if os.path.exists(os.path.join(session_dir, 'beh_df.csv')) and overwrite==False:

        beh_df = pd.read_csv(os.path.join(session_dir, 'beh_df.csv'))

    else:
        beh_info = load_behavior(beh_path)
        beh_df = make_beh_df(beh_info)

    return beh_df, spks_info

### ---fx called by load_and_wrangle---
def load_behavior(beh_path):
    """
    This function loads the behavior data from the protocol_info.mat file in the
    directory that is passed into it

    inputs
    ------
    beh_path  : string, path to .mat file with behavior info

    returns
    ------
    beh_info  : ndarray, with behavior .mat structure extracted
    """
    beh_info = load_nested_mat(beh_path) # this function is custom to deal with nested structures
    beh_info = beh_info["behS"]

    return beh_info


def load_spks(spks_path):
    """
    This function loads the behavior data from the ksphy_cluster_info_foranalys.mat
    file in the directory that is passed into it

    inputs
    ------
    spks_path : string, path to .mat file with spiking info extracted from kilosort/raw data

    returns
    ------
    spks_info : ndarray, with spk info from .mat file
    """
    spks_info = spio.loadmat(spks_path)
    spks_info = spks_info['PWMspkS'][0]

    return spks_info


def make_beh_df(beh_info):

    """
    Data wrangling function to take dictionary from load_session_info() and putting it into a
    tidy dataframe for analysis

    inputs
    ------
    beh_info : ndarray, extracted .mat structure from load_session_info()

    returns
    -------
    beh_df   : df (ntrials x items), tidy data frame with behavior information & some relabeling
    """

    # initialize data frame
    beh_df = pd.DataFrame()

    # assign trail n values
    beh_df['trial_num'] = np.arange(1, beh_info['n_completed_trials'] + 1)

    # rename previous side to be L/R
    prev_side_adj = np.roll(beh_info['prev_side'],1) # n-1 trial info
    prev_side_adj = np.where(prev_side_adj == 114, 'RIGHT', 'LEFT' )
    prev_side_adj[0] = 'NaN' # trial 0 doesn't have a previous

    # turn hit info to strings
    beh_df['hit_hist'] = beh_info['hit_history']
    beh_df['hit_hist'] = beh_df['hit_hist'].mask(beh_df['hit_hist'] == 1.0, "hit")
    beh_df['hit_hist'] = beh_df['hit_hist'].mask(beh_df['hit_hist'] == 0.0, "miss")
    beh_df['hit_hist'][beh_df['hit_hist'].isnull()] = "viol"

    # get n_trial length items into df
    beh_df['delay'] = beh_info['delay']
    beh_df['pair_hist'] = beh_info['pair_history']
    beh_df['correct_side'] = beh_info['correct_side']
    beh_df['prev_side'] = prev_side_adj
    beh_df['aud1_sigma'] = beh_info['aud1_sigma']
    beh_df['aud2_sigma'] = beh_info['aud2_sigma']

    # extract parsed events/state machine info for each trial
    parsed_events_dict = beh_info["parsed_events"]

    # initilize space
    c_poke = np.zeros((len(parsed_events_dict)))
    hit_state = np.zeros((len(parsed_events_dict)))
    aud1_on = np.zeros((len(parsed_events_dict)))
    aud1_off = np.zeros((len(parsed_events_dict)))
    aud2_on = np.zeros((len(parsed_events_dict)))
    aud2_off = np.zeros((len(parsed_events_dict)))
    end_state = np.zeros((len(parsed_events_dict)))

    # iterate over items from state matrix
    for trial in range(len(parsed_events_dict)):

        # every trial has a center poke & end_state
        c_poke[trial] = parsed_events_dict[trial]['states']['cp'][0]
        end_state[trial] = parsed_events_dict[trial]['states']['check_next_trial_ready'][1]

        # not all trials will have sound/hit time/etc, pull out info for non-violated
        if beh_df['hit_hist'][trial] == 'viol':
            hit_state[trial] = float("NaN")
            aud1_on[trial]   = float("NaN")
            aud1_off[trial]  = float("NaN")
            aud2_on[trial]   = float("NaN")
            aud2_off[trial]  = float("NaN")

        elif beh_df['hit_hist'][trial] == 'hit':
            hit_state[trial] = parsed_events_dict[trial]['states']['hit_state'][0]
            aud1_on[trial]   = parsed_events_dict[trial]['waves']['stimAUD1'][0]
            aud1_off[trial]  = parsed_events_dict[trial]['waves']['stimAUD1'][1]
            aud2_on[trial]   = parsed_events_dict[trial]['waves']['stimAUD2'][0]
            aud2_off[trial]  = parsed_events_dict[trial]['waves']['stimAUD2'][1]

        elif beh_df['hit_hist'][trial] == 'miss':
            hit_state[trial] = parsed_events_dict[trial]['states']['second_hit_state'][0]
            aud1_on[trial]   = parsed_events_dict[trial]['waves']['stimAUD1'][0]
            aud1_off[trial]  = parsed_events_dict[trial]['waves']['stimAUD1'][1]
            aud2_on[trial]   = parsed_events_dict[trial]['waves']['stimAUD2'][0]
            aud2_off[trial]  = parsed_events_dict[trial]['waves']['stimAUD2'][1]

        else:
            raise Exception('hit_hist doesn''t appear to have correct structure (hit, miss, viol)')


    # add to df
    beh_df['c_poke']    = c_poke
    beh_df['end_state'] = end_state
    beh_df['hit_state'] = hit_state
    beh_df['aud1_on']   = aud1_on
    beh_df['aud1_off']  = aud1_off
    beh_df['aud2_on']   = aud2_on
    beh_df['aud2_off']  = aud2_off
    beh_df['end_state'] = end_state

    return beh_df
### --- end of fx called by load_and_wrangle
