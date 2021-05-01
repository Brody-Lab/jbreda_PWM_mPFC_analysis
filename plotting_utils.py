
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
from scipy import stats
# stored one repo up in my fork of Spykes
from spykes.spykes.plot.neurovis import NeuroVis


"PSTHs- Gaussain"

def PSTH_gaussain(event_aligned_spks, event_aligned_windows, event, df,
                  conditions=None, bin_size=0.001, mu=0, sigma=0.150):

    """
    Function for computing PSTH information w/ guassain smoothing

    params:
    -------
    event_aligned_spks    : dict, event_aligned spike times for a single neuron & event from align_neuron_to_events()
    event_aligned_windows : dict, windows of time used in align_neuron_to_events()
    event                 : str, event to align to (this is a key for the two alignment dictionaries)
    df                    : df, behavior information (usually filtered) to use for determining condition trials
    conditions            : str, defualt = None, name of column in df to split by
    bin_size              : int, default = 0.001, time in s used to binarize spike trian
    mu                    : int, default = 0, mean of guassian kernal in seconds
    sigma                 : int, faults = 0.150, std of guassain kernal in seconds

    returns
    -------
    psth                  : dict, containing smoothed data for a single neuron and event to be used to use in
                            plotting or further analysis

    """

    ## biniarize spike train
    binarized_spks = binarize_event(event_aligned_spks[event], event_aligned_windows[event],
                                    bin_size=bin_size)

    ## make kernal
    x = np.linspace(-5,5, 350)
    kernal = make_gaussian_kernal(x, mu=mu, sigma=sigma)

    ## get a set of binary indicators to split by conditions (if present)
    trials = dict()
    if conditions:
        for cond_id in np.sort(df[conditions].unique()):

            trials[cond_id] = np.where((df[conditions] == cond_id).values)[0]
    else:
        trials['all'] = np.where(np.ones(len(df)))[0]

    ## iterate over each condition, smooth with kernal and save to psth dict
    psth = {'event' : event, 'conditions' : conditions, 'n': [], 'time' : [], 'data' : {}, 'mean' : {}, 'sem' :{}}

    for cond_id, idxs in trials.items():

        # grab the trials for each condition
        selected_trials = binarized_spks[idxs]

        # convolve
        mean, sem, data = smooth_trials(selected_trials, kernal, summary=True)

        # append
        psth['n'].append(len(selected_trials))
        psth['data'][cond_id] = data
        psth['mean'][cond_id] = mean
        psth['sem'][cond_id] = sem
        psth['time'].append(np.linspace(event_aligned_windows[event][0],
                                    event_aligned_windows[event][1],
                                    len(psth['mean'][cond_id])))

    return psth

def binarize_event(event_aligned_spks, window, bin_size):

    """
    Function for taking event centered spike times and binarizing them for PSTHs

    params
    ------
    event_aligned_spks : event_aligned spike times for a single neuron & event from align_neuron_to_events()
    window : window of time around event to binarize (NOTE: this is limited by the windows in simuli_align()
    bin_size : int, binsize in s to use for determining if there is a spike or not

    returns
    -------
    binarized_spks : list, (1 x n_trials) with binarized spike counts for a single event and cell

    notes
    -----
    after checking a few neruons, avereage ISI is 10-50 ms, bin_size should be < 10 and > 1 ms to balance
    accuracy and speed
    """

    # initialize
    binarized_spks = []
    half_bin = bin_size / 2
    bin_centers = np.arange((window[0] * 0.001) + half_bin, (window[1] * 0.001), bin_size)
    n_bins = len(bin_centers)

    # iterate over each trial for an event
    for itrial in range(len(event_aligned_spks)):

        binarized_trial = np.zeros((n_bins))

        # iterate over each time bin in a trial
        for ibin in range(n_bins):

            # for each of the spikes in the ith trial, do any fit in the ith bin?
            spk_in_ibin = np.logical_and(event_aligned_spks[itrial] >= (bin_centers[ibin] - half_bin),
                                        event_aligned_spks[itrial] <= (bin_centers[ibin] + half_bin))

            # if yes, report it
            if np.sum(spk_in_ibin) > 0:
                binarized_trial[ibin] = 1

        binarized_spks.append(binarized_trial)

    return np.array(binarized_spks)

def gaussian(x, mu, sigma):
    "Quick fx for guassian distribution"
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1/2 * ((x - mu)/sigma)**2)

def make_gaussian_kernal(x, mu, sigma):
    "Qucik function for making a guassian kernal wtih specified mean and std dev"

    kernal = gaussian(x, mu, sigma)
    kernal_normalized = kernal/np.sum(kernal) # area = 1

    return kernal_normalized

def smooth_trial(binarized_trial, kernal):
    """
    Function that convolved spikes from a single trial with a kernal. Because the sigma is large,
    it will drop any signal that is 0 and replace with nan's because this is where masking occured
    !!NOTE!! this should be updated for new animals after W122
    """
    # multply by 1000 to get to spks/second (as opposed to spks/ms)
    smoothed = np.convolve(binarized_trial, kernal, mode = 'same') * 1000
    smoothed_remove_masking = np.where(smoothed == 0, np.nan, smoothed)
    return smoothed_remove_masking

def smooth_trials(binarized_trials, kernal, summary):

    """
    Function for smoothing event cented trials (for the whole session, or a condiiton)
    given a kernal.

    params
    ------
    binarized_trials : list, output from binarize_event()
    kernal : arr, kernal to convole with
    summary : bool, whether to calculate and return summary information for trials

    returns
    -------
    smoothed_mean : arr, mean of smooothed trials
    smoothed_sem : arr, standard error of mean of smoothed trials
    smoothed_trial : arr, smoothed (aka convolved) output for each trial
    """

    smoothed_trials = []

    for trial in binarized_trials:
        smoothed_trials.append(smooth_trial(trial, kernal))

    if summary:
        smoothed_mean = np.nanmean(np.array(smoothed_trials), axis=0)
        smoothed_sem = stats.sem(np.array(smoothed_trials), axis=0, nan_policy='omit')

        return smoothed_mean, smoothed_sem, smoothed_trials

    else:
        return np.array(smoothed_trials)

"PSTHs- Boxcar"

def PSTH_boxcar(event_aligned_spks, event_aligned_windows, event, df,
                  conditions=None, bin_size=0.100, masking=True):

    """
    Function for computing PSTH information w/ guassain smoothing

    params:
    -------
    event_aligned_spks    : dict, event_aligned spike times for a single neuron & event from align_neuron_to_events()
    event_aligned_windows : dict, windows of time used in align_neuron_to_events()
    event                 : str, event to align to (this is a key for the two alignment dictionaries)
    df                    : df, behavior information (usually filtered) to use for determining condition trials
    conditions            : str, defualt = None, name of column in df to split by
    bin_size              : int, default = 0.100, time in s to make box car
    masking               : bool, default = True, if 0s should be set to nans due to pre-kilosort masking

    returns
    -------
    psth                  : dict, containing summed data for a single neuron and event to be used to use in
                            plotting or further analysis

    """

    assert len(df) == len(event_aligned_spks[event]), "Trial and spike lists of different length"

    ## get a set of binary indicators to split by conditions (if present)
    trials = dict()

    if conditions:
        for cond_id in np.sort(df[conditions].unique()):

            trials[cond_id] = np.where((df[conditions] == cond_id).values)[0]
    else:
        trials['all'] = trials['all'] = np.where(np.ones(len(df)))[0]

    ## iterate over each condition, count, summarize and save to psth dict
    psth = {'event' : event, 'conditions' : conditions, 'n': [], 'time' : [], 'data' : {}, 'mean' : {}, 'sem' :{}}

    for cond_id, idxs in trials.items():

        print(cond_id, idxs)

        # grab the trials for each condition
        selected_trials = np.array(event_aligned_spks[event], dtype=object)[idxs]
        # print(selected_trials)


        # count
        counted_spks = get_spike_counts(selected_trials, event_aligned_windows[event], bin_size=bin_size)
        print('counted spks:', counted_spks)
        # summarize
        counted_spks_masked, mean, sem = summarize_spike_counts(counted_spks, masking=masking)

        # append (divide by bin_size to get in Hz)
        psth['n'].append(len(selected_trials))
        psth['data'][cond_id] = counted_spks_masked /  bin_size
        psth['mean'][cond_id] = mean / bin_size
        psth['sem'][cond_id] = sem / bin_size
        psth['time'].append(np.linspace(event_aligned_windows[event][0],
                                    event_aligned_windows[event][1],
                                    len(psth['mean'][cond_id])))

    return psth

def get_spike_counts(event_aligned_spks, window, bin_size):

    """
    Function for taking event centered spike times and counting them using sliding boxcar

    params
    ------
    event_aligned_spks : event_aligned spike times for a single neuron & event from align_neuron_to_events()
    window : window of time around event to analyze (NOTE: this is limited by the windows in simuli_align()
    bin_size : int, binsize in s to use for boxcar

    returns
    -------
    counted_spks : list, (1 x n_trials) with spike counts for a single event and cell given bin_size

    """
    counted_spks = []
    half_bin = bin_size / 2
    bin_centers = np.arange((window[0] * 0.001) + half_bin, (window[1] * 0.001), bin_size)
    n_bins = len(bin_centers)

    # iterate over each trial for an event
    for itrial in range(len(event_aligned_spks)):

        counted_trial = np.zeros((n_bins))

         # iterate over each time bin in a trial
        for ibin in range(n_bins):

            # for each of the spikes in the ith trial, how many fit into the ith bin?
            n_spikes_in_bin = np.logical_and(event_aligned_spks[itrial] >= (bin_centers[ibin] - half_bin),
                                            event_aligned_spks[itrial] <= (bin_centers[ibin] + half_bin))

            counted_trial[ibin] = np.sum(n_spikes_in_bin)

        counted_spks.append(counted_trial)

    return np.array(counted_spks)

def summarize_spike_counts(counted_spks, masking=True):
    "Quick fx for getting summary information and nan-packing of spike count data"

    if masking:

        # replace 0s with nans
        counted_spks = np.where(counted_spks == 0, np.nan, counted_spks)

    counted_mean = np.nanmean(counted_spks, axis = 0)
    counted_sem = stats.sem(counted_spks, axis=0, nan_policy='omit')

    return counted_spks, counted_mean, counted_sem

""" ITEMS BELOW USED FOR SPYKES """


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
