
#!/usr/bin/env python

"""
Utils to be used for creating objects in for plotting with Spykes and iterative
plotting w/ modified Spykes functions.
"""

# libraries
import sys; sys.path.insert(0, '..') # if you don't find it here, look one above
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.io as spio
# stored one repo up in my fork of Spykes
from spykes.spykes.plot.neurovis import NeuroVis

def split_spks_by_trial(spk_times, df):
    """
    Function that converts spikes for a whole trial (1d) into array where rows are trials
    and columns are spk times within a trial for a neuron, added to dict with other trial level info

    inputs
    ------
    spk_times : arr-like, 1d containing spk times in seconds for a single neuron for a session
    df        : df, containing behavior information for a trial

    returns
    -------
    trial_spks_dict : dict, containing trial number (w/ respect to whole session), start and stop times
                      for a trial and the spks occuring in a a trial in s
    """

    # initialize
    trial_spks_dict = {}
    trial_nums = []
    trial_starts = []
    trial_stops = []
    trial_spks = []

    # make counter because dataframe idx is not continous
    # instead, it represents overall trial number with respect to start of session
    itrial = 0

    for idx, trial in df.iterrows():

        # grab start time 500 ms before c poke
        trial_start = trial['c_poke'] -0.5

        # create end time so delays have same length
        if trial['delay'] == 2:
            trial_end = trial_start + 5.5
        elif trial['delay'] == 4:
            trial_end = trial_start + 7.7
        elif trial['delay'] == 6:
            trial_end = trial_start + 9.9
        else:
            print('Trial delay is not 2, 4 or 6 seconds- this is not correct')

        # append trial info
        trial_nums.append(idx)
        trial_starts.append(trial_start)
        trial_stops.append(trial_end)

        # grab spike time within trial window & append them
        trial_spks.append(spk_times[np.logical_and(spk_times >= trial_start, spk_times <= trial_end)])


    trial_spks_dict['trial_nums'] = trial_nums
    trial_spks_dict['trial_times'] = trial_starts, trial_stops
    trial_spks_dict['trial_spks'] = trial_spks

    return trial_spks_dict

def get_spike_counts(trial_spks_dict, bin_size, mode, trial_len='same'):

    """
    Function that takes spike time information for a session split by trials and binarizes or counts the
    spike information for further analyses

    Inputs
    ------
    trial_spks_dict     : dict, created by split_spks_by_trial()
    bin_size            : int, size of sliding bin to use in seconds
    mode                : str, whether binirze spike counts for a bin ('binary)', or count them ('count')
    trial_len           : str, whether to pack trials of different lenghts with 0s to make all the same
                          length ('same') or use the ending time ('end'), default='same'

    Returns
    -------
    session_spk_binary  : array_like, N trials X trial_len, biniarize spike train for each trial

    Notes
    -----
    For binarizing data, a 0.001 s window does quite well, but there will still be cases of 2-3 spikes in a bin
    only being counted as 1 spike. This appears to be quite rare on the scale of 2-5 bins on 25% of trials
    """

    # initialize & extract
    session_spk_counts = []
    spks = trial_spks_dict['trial_spks']

    # iterate over spikes for each trial
    for itrial, trial_spks in enumerate(tqdm(spks)):

        # grab start time
        t_start = trial_spks_dict['trial_times'][0][itrial]

        # determine which end time to use given inputs
        # NOTE this assumes 6s delays are included so maximum t_len = 9.9 seconds
        if trial_len == 'same':
            t_end = t_start + 9.9
        elif trial_len == 'end':
            t_end = trial_spks_dict['trial_times'][1][itrial]
        else:
            print('This is not a valid_trial length type')

        # initialize bin structure to include all bins form [0 to bin_size]
        half_bin = bin_size / 2
        bin_centers = np.arange(t_start + half_bin, t_end, bin_size)
        n_bins = len(bin_centers)

        # updated for each trial
        trial_spk_counts = np.zeros((n_bins))

        for ibin in range(n_bins):
            # for each of the spikes in the ith trial, do any fit in the ith bin?
            if mode == 'binary':
                spike_in_ibin = np.logical_and(spks[itrial] >= (bin_centers[ibin] - half_bin),
                                               spks[itrial] <= (bin_centers[ibin] + half_bin))

                # if there is a spike in the bin, report it
                if np.sum(spike_in_ibin) > 0:
                    trial_spk_counts[ibin] = 1

            elif mode == 'count':
                # for each of the spikes in the ith trial, how many fit into the ith bin?
                n_spike_in_bin = np.sum(np.logical_and(spks[itrial] >= (bin_centers[ibin] - half_bin),
                                                       spks[itrial] <= (bin_centers[ibin] + half_bin)))

                trial_spk_counts[ibin] = n_spike_in_bin

        session_spk_counts.append(trial_spk_counts)

    return session_spk_counts


def initiate_neurons(spks_dict):
    """
    This function takes spks_dict along with session data and unpacks spike times
    in finite state machine/behavior time & turns into NeuroVis object.

    inputs
    ------
    spks_dict : dict, with spikes .mat structures extracted
    sess_data : str, used for naming the objects

    returns
    -------
    neuron_list : list, with each item being a NeuroVis object pertaining to a
                  neuron in the session
    """

    spk_in_fsm_time = spks_dict["spk_times"] # fsm = behavior time
    sess_date = spks_dict['date']
    neuron_list = []

    for neuron in range(len(spk_in_fsm_time)):
        spk_times = spk_in_fsm_time[neuron]

        # instantiate neuron
        neuron = NeuroVis(spk_times, name = '{} {}'.format(neuron, sess_date))
        neuron_list.append(neuron)

    return neuron_list

"get_raster for all neurons and assigned events"


def get_neuron_rasters(neurons, events, windows, bndl_dfs, df_names, conditions=None, binsize = 50):
    """
    This function can be used to get rasters for multiple neurons across multiple events
    with specific time windows for each event

    inputs:
    -------
    neurons   : NeuroVis object, N neurons long
    events    : list, event name in strings from your behavior df you want to align to
    windows   : list, time window in ms to grab for each event
    bndl_dfs  : dict, containing df for each bndl that has a mask created by
                make_unmasked_dfs()
    df_names  : list, Ncells long containing dict keys to access df for each
                cell  created by make_unmasked_dfs()
    returns:
    -------
    neuron_rasters : list, raster dictionaries stored by [neuron][event] for plotting"""

    #initialize storage ([Neuron][Event])
    neuron_rasters = []

    # iterate over each neuron and event
    for neuron in range(len(neurons)):
        rasters = []

        for event in range(len(events)):

        # create raster dictionary
            raster = neurons[neuron].get_raster(event = events[event], conditions=conditions,
                                                df = bndl_dfs[df_names[neuron]],
                                                window=windows[event], plot=False,
                                                binsize=binsize)
            rasters.append(raster)

        neuron_rasters.append(rasters)

    return neuron_rasters


def get_neuron_psths(neurons, events, windows, bndl_dfs, df_names, conditions=None, binsize=50):
    """
    This function can be used to get psths for multiple neurons across multiple events
    with specific time windows for each event

    inputs:
    -------
    neurons    : NeuroVis object, N neurons long
    events     : list, event name in strings from behavior df you want to align to
    conditions : str, condition namefrom behavior df to split by (e.g. hit Y/N)
    windows    : list, time window in ms to grab for each event
    binsize    : int, binsize in ms
    bndl_dfs  : dict, containing df for each bndl that has a mask created by
                make_unmasked_dfs()
    df_names  : list, Ncells long containing dict keys to access df for each
                cell  created by make_unmasked_dfs()

    returns:
    -------
    neuron_psths : list, psths dicts stored by [neuron][event] for plotting"""


    #initialize storage ([Neuron][Event])
    neuron_psths = []

    # iterate over each neuron and event
    for neuron in range(len(neurons)):
        psths = []

        for event in range(len(events)):

        # create psth dictionary
            psth = neurons[neuron].get_psth(event=events[event], df=bndl_dfs[df_names[neuron]],
                          window=windows[event], conditions=conditions, binsize=binsize, plot=False,
                          event_name=events[event])

            psths.append(psth)

        neuron_psths.append(psths)

    return neuron_psths
