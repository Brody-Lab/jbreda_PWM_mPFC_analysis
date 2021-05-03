
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
from scipy.ndimage import gaussian_filter1d
import statsmodels.api as sm
# stored one repo up in my fork of Spykes
from spykes.spykes.plot.neurovis import NeuroVis

delay_colors =['#e1e1e1','#e1e1e1', '#929292', '#4b4b4b', '#1e1e1e']


"PSTHs- Gaussain"

def PSTH_gaussain(event_aligned_spks, event_aligned_windows, event, df, conditions=None,
                  bin_size=0.001, sigma=150):

    """
    Function for computing PSTH information w/ guassain smoothing

    params:
    -------
    event_aligned_spks    : dict, event_aligned spike times for a single neuron & event from align_neuron_to_events()
    event_aligned_windows : dict, windows of time used in align_neuron_to_events()
    event                 : str, event to align to (this is a key for the two alignment dictionaries)
    df                    : df, behavior information (usually filtered) to use for determining condition trials
    conditions            : str, optional- defualt = None, name of column in df to split by
    bin_size              : int, optional- default = 0.001, time in s used to binarize spike trian
    x                     : arr, optional- values used to make the kernal w/ mu and sigma
    mu                    : int, optional- default = 0, mean of guassian kernal in seconds
    sigma                 : int, optional- default = 0.150, std of guassain kernal in seconds

    returns
    -------
    psth                  : dict, containing smoothed data for a single neuron and event to be used to use in
                            plotting or further analysis

    """

    ## biniarize spike train
    binarized_spks = binarize_event(event_aligned_spks[event], event_aligned_windows[event],
                                    bin_size=bin_size)


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
        mean, sem, data = smooth_trials(selected_trials, sigma=sigma, summary=True)

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

def smooth_trial(binarized_trial, sigma):
    """
    Function that convolved spikes from a single trial with a kernal. Because the sigma is large,
    it will drop any signal that is 0 and replace with nan's because this is where masking occured
    !!NOTE!! this should be updated for new animals after W122

    sigma : std deviation of gaussian in ms
    """
    # multply by 1000 to get to spks/second (as opposed to spks/ms)
    smoothed = gaussian_filter1d(binarized_trial, sigma, mode='wrap') * 1000
    smoothed_remove_masking = np.where(smoothed == 0, np.nan, smoothed)

    return smoothed_remove_masking

def smooth_trials(binarized_trials, sigma, summary):

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
        smoothed_trials.append(smooth_trial(trial, sigma))

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

        # grab the trials for each condition
        selected_trials = np.array(event_aligned_spks[event], dtype=object)[idxs]

        # count
        counted_spks = get_spike_counts(selected_trials, event_aligned_windows[event], bin_size=bin_size)

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

def plot_psth(psth, axis=None, title=None, xlim=None, ylim=None,
              legend=False, stimulus_bar=None):

    """
    Function for plotting PSTHs

    params:
    -------
    psth  : dict, output from PSTH_boxcar() or PSTH_gaussain()
    axis  : optional- default = gca(), which axis to plot to
    title : str, optional- default = None, title of plot
    xlim  : arr, optional- default = None, x limits of plot
    ylim  : arr, optional- default = None, y limits of plot
    legens : bool, optional- default = False, to turn legend on
    stimulus_bar : str, optional- default = None, 'sound on' or 'sound off'

    returns:
    --------
    none, only plots
    """

    # set axes
    if axis:
        ax = axis
    else:
        ax = plt.gca()

    # pull of of dictionary for ease
    mean = psth['mean']
    sem = psth['sem']
    time = psth['time'][0] # all times are the same, just grab one

    # plot by condition
    for idx, cond_id in enumerate(mean.keys()):

        ax.plot(time, mean[cond_id], color=delay_colors[idx])


        ax.fill_between(time, mean[cond_id] - sem[cond_id],
                       mean[cond_id] + sem[cond_id], alpha = 0.2,
                       color=delay_colors[idx], label=cond_id)

    ax.axvspan(0,0.0001, color = 'black')

    # aesthetics
    ax.set(xlabel = 'Time (ms)', ylabel = 'Firing Rate (Hz)')

    if title:
        ax.set_title(title)

    if xlim:
        ax.set_xlim(xlim)

    if ylim:
        ax.set_ylim(ylim)

    else:
        scale = 0.01
        y_min = (1 - scale) * np.nanmin([np.min(mean[cond_id]) for cond_id in mean.keys()])
        y_max = (1 - scale) * np.nanmax([np.max(mean[cond_id]) for cond_id in mean.keys()])

    if legend:
        ax.legend(frameon=False, bbox_to_anchor=(1, 1))

    if stimulus_bar == 'sound on':
        ax.axvspan(0, 400, alpha=0.2, color='grey')

    elif stimulus_bar == 'sound off':
        ax.axvspan(-400, 0, alpha =0.2, color='grey')

    sns.despine()

def fr_by_loudness_df(psth, neuron_id):

    """
    Create df with average firing rate information by first sound loudness for
     further plotting/analysis

    params:
    -------
    psth      : dict, output from psth_gaussain() or psth_boxcar() with firing rate
                information
    neuron_id : str, session data and neuron idx for labeling

    returns
    -------
    df        : df, n trials long with first sound loduness and average firing
                rate during delay period"""

    conds = []
    mean_fr_by_cond = []

    for key in psth['data'].keys():
        for trial_psth in psth['data'][key]:

            # update loudness values
            conds.append(float(key.replace('*','')))

            # get mean for each trial during only the delay period
            mean_fr_by_cond.append(np.nanmean(trial_psth[150:-150]))

    ids = [neuron_id] * len(conds)

    df = pd.DataFrame({'firing_rate' : mean_fr_by_cond, 'condition' : conds, 'neuron_id' : ids})

    return df

def simple_regplot( x, y, n_std=2, n_pts=100, ax=None, scatter_kws=None, line_kws=None,
    ci_kws=None, title=None, xlabel=None, ylabel=None):

    """ Draw a regression line with error interval & save output from stats model.
    Modified from :https://stackoverflow.com/questions/22852244/how-to-get-the-
    numerical-fitting-results-when-plotting-a-regression-in-seaborn """

    ax = plt.gca() if ax is None else ax

    # calculate best-fit line and interval
    x_fit = sm.add_constant(x)
    fit_results = sm.OLS(y, x_fit, missing='drop').fit()

    eval_x = sm.add_constant(np.linspace(np.min(x), np.max(x), n_pts))
    pred = fit_results.get_prediction(eval_x)

    # draw the fit line and error interval
    ci_kws = {} if ci_kws is None else ci_kws
    ax.fill_between(
        eval_x[:, 1],
        pred.predicted_mean - n_std * pred.se_mean,
        pred.predicted_mean + n_std * pred.se_mean,
        alpha=0.5,
        **ci_kws,
    )
    line_kws = {} if line_kws is None else line_kws
    h = ax.plot(eval_x[:, 1], pred.predicted_mean, **line_kws)

    # draw the scatterplot
    scatter_kws = {} if scatter_kws is None else scatter_kws
    ax.scatter(x, y, color=h[0].get_color(), **scatter_kws)

    # add information
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    sns.despine()

    return fit_results

def loudness_regression(df, ax=None):

    """
    Simple linear regression on loudness by firing rate with visualization

    params
    ------
    df : df, dataframe for a neuron created by fr_by_loudness_df()

    returns
    -------
    stats_df : df, r^2 and p-value output for fr ~ condition regression

    plots
    -----
    firing rate by condition performed on all trials, mean denoted in black
    """

    neuron_id = df['neuron_id'][0]
    ax = plt.gca() if ax is None else ax


    fit = simple_regplot(df['condition'], df['firing_rate'], ax=ax, title=neuron_id,
                         ci_kws={'color':'grey'}, line_kws={'color':'grey'},
                         scatter_kws={'alpha':0.5}, xlabel='Loudness (dB)',
                         ylabel='Firing Rate (Hz)')

    # get summary info & mark black on plot for visual
    df_summary = df.groupby(['condition'], as_index=False).mean()
    ax.scatter(df_summary['condition'], df_summary['firing_rate'], color='black')

    # add in statistics & save them out
    ax.text(0.055, df['firing_rate'].max(),f'$R^2$ = {fit.rsquared:0.2f} \n p = {fit.pvalues[1]:0.2f}')

    stats_df = pd.DataFrame({'neuron_id' : neuron_id,
                            'rsqaured' : fit.rsquared,
                            'pvalue' : fit.pvalues[1]}, index=[0])
    return stats_df
