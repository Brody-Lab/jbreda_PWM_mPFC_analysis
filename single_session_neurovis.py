 #!/usr/bin/env python
# coding: utf-8
"""
Written by JRB 2021-02-22

Purpose
-------
This script takes a session name as an argument and for each neuron in that session
and event specififed below, makes a figure w/ raster, psth and average waveform.

Notes
-----
Currently, this is written for plotting 2 second delay hit trials
"""
# libraries
import sys; sys.path.insert(0, '..') # if you don't find it here, look one above
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.io as spio
import pickle
# stored one repo up in my fork of Spykes
from spykes.spykes.plot.neurovis import NeuroVis
from io_utils import *
from plotting_utils import *
sns.set_context("talk")

# --- IO ---
# base paths/names
sess_name=sys.argv[1]
base_path  = 'Y:\jbreda\ephys\post_sort_analysis\sorted_pre_bdata'
beh_mat   = 'protocol_info.mat'
spks_mat  = 'ksphy_clusters_foranalysis.mat'

# create paths
sess_path = os.path.join(base_path, sess_name)
beh_path  = os.path.join(sess_path, beh_mat)
spks_path = os.path.join(sess_path, spks_mat)

fig_save_path = os.path.join(os.getcwd(), 'figures', 'neurovis', sess_name)
if os.path.exists(fig_save_path):
    print('fig save path exists, will not create')
else:
    os.mkdir(os.path.join(fig_save_path))

# load & wrangle
beh_df, spks_dict = load_and_wrangle(beh_path, spks_path, overwrite=True)


# --- Initialize for plotting ---
# filter dataframe for 2 second hit trials
beh_df_d2_h = beh_df[beh_df['hit_hist'] == 'hit']

# deal with masking
bndl_dfs, df_names = deal_with_masking(spks_dict, beh_df_d2_h, sess_path)

# assign events, windows and conditions
events = ['c_poke', 'aud1_on', 'aud2_on']
windows = [[-300,700], [-500,1000], [-500, 500]]
#[-1000,6300]]
condition = 'delay'

# get neurons into NeuroVis objects
neurons = initiate_neurons(spks_dict)

# create
neuron_rasters = get_neuron_rasters(neurons, events, windows, bndl_dfs, df_names)
neuron_psths = get_neuron_psths(neurons, events, windows, bndl_dfs, df_names, conditions=None)


# --- Plot ---
# iterate over neurons
print('Plotting Start')
for nn in range(len(neurons)):

    # grab mean waveform for the neuron
    mean_wav = spks_dict['mean_wav'][nn]

    # iterate over event
    for ee in range(len(events)):

        # initilaize plot
        fig = plt.figure(figsize=(18,12))
        ax1 = plt.subplot2grid((4, 3), (0, 0), rowspan=2, colspan=2)
        ax2 = plt.subplot2grid((4, 3), (2, 0), rowspan=2, colspan=2)
        ax3 = plt.subplot2grid((4, 3), (1, 2), rowspan=2, colspan=1)

        # plot raster and psth for each neuron, event
        neurons[nn].plot_raster(neuron_rasters[nn][ee], axis = ax1, event_name=events[ee], cmap="Greys")
        neurons[nn].plot_psth(neuron_psths[nn][ee],axis=ax2, event_name=events[ee])

        # plot average waveform from active tetrode
        for tt in range(4):
            ax3.plot(mean_wav[tt])
            ax3.set_title('Average Waveform')

        fig_name = "{sess}_neuron_{neuron}_{event}_{ee}".format(sess = sess_name,
                                                                neuron = nn, event = events[ee],
                                                               ee = ee)
        # format & save out
        plt.tight_layout()
        plt.savefig(os.path.join(fig_save_path, fig_name))
        plt.close("all")
print('Plotting End')
