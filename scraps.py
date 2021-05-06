def gaussian(x, mu, sigma):
    "Quick fx for guassian distribution"
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1/2 * ((x - mu)/sigma)**2)

def make_gaussian_kernal(x, mu, sigma):
    "Qucik function for making a guassian kernal wtih specified mean and std dev"

    kernal = gaussian(x, mu, sigma)
    kernal_normalized = kernal/np.sum(kernal) # area = 1


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


def run_delay_analysis_for_multiple_events(sess_name, sess_aligned, align_windows, events, dfs,
                       fig_save_path, sess_path):
    """
    Wrapper function for analyze_and_plot_loudness() to allow for analysis on multiple events
    in a session with varying trials (ie, 2 and 4s delay)

    Params
    ------
    sess_name : str, name of the session for id in plotting
    sess_aligned : dict, nested with dicts for each neuron giving alignment times for events
                   in the trial output from event_align_session()
    sess_windows : dict, with timing information from event_align_session(), for a single neuron
                  (because it's the same for all neurons in a session)
    events        : list, events in which you want to align to, used as a key for sess_algined,
                   sess_windows. E.g. ['delay2s', 'delay4s']
    dfs           : list, behavior data frames to use based on your alignment events. For example,
                    if aligning to `delay2s`, your df should contain only 2s trials
    fig_save_path : str, where you want to save out the psth and loudness regression figures

    sess_path     : str, path to session data where analysis .csv files will be stored with
                    regression output

    plots & saves
    -----
    - psth for each neuron for each event, split by the loduness of the first sound
    - OLS regression of first sound on firing rate for each neuron
    - table with average firing rate for each trial used in regression for each neuron
    - table with p-value and r^2 for each neuron regression

    """

    for event, df in zip(events, dfs):

        analyze_and_plot_loudness(sess_name, sess_aligned, align_windows,
                                  event, df, fig_save_path, sess_path)





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
