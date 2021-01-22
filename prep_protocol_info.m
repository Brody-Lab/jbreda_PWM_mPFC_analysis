%% Prep protocol info
% The purpose of this script is to use all the information extracted in
% find_wireless_sess_J to get into bdata & extract protcol information for
% a single rat/session to be used later on for ephys analysis in python.
% ------------------------------
%% Initialize Inputs

sess_name = 'data_sdc_20190902_145404_fromSD';
sorted_dir =  'Y:\jbreda\ephys\post_sort_analysis\sorted_pre_bdata';
sess_dir = fullfile(sorted_dir, sess_name);
save_name  = fullfile(sess_dir,'protocol_info.mat')

%% Load up
% load spkS & alignmnet info
[ ~ , PWMspkS] = get_ksphy_results_J(sess_name, 'overwrite', false);
 
% load behavior data & dohble check we are in the right place
beh_dat = load(PWMspkS(1).behav_session.goodpath);
assert(beh_dat.saved.PWM_sessid == PWMspkS(1).sessid)

% single items
behS.ratname            = beh_dat.saved.SavingSection_ratname;
behS.ratmass            = beh_dat.saved.AdLibGUISection_masses;
behS.sessid             = beh_dat.saved.PWM_sessid;
behS.n_started_trials   = beh_dat.saved.ProtocolsSection_n_started_trials;
behS.n_completed_trials = beh_dat.saved.ProtocolsSection_n_completed_trials;
behS.hit_frac           = beh_dat.saved.RewardsSection_hit_frac;
behS.Lhit_frac          = beh_dat.saved.RewardsSection_Left_hit_frac;
behS.Rhit_frac          = beh_dat.saved.RewardsSection_Right_hit_frac;

% psychometrics
n_sound_pairs = 14;
for n_pair = 1:n_sound_pairs
    base = 'PWMSection_';
    prob = sprintf('probClass%i', n_pair);
    perf = sprintf('perfClass%i', n_pair);
    ntrials = sprintf('nTrialsClass%i', n_pair);
   
    behS.psycho(n_pair).prob    = beh_dat.saved.(strcat(base, prob));
    behS.psycho(n_pair).perf    = beh_dat.saved.(strcat(base, perf));
    behS.psycho(n_pair).ntrials = beh_dat.saved.(strcat(base, ntrials));
end

behS.soundpairs                 = beh_dat.saved.PWMSection_thesepairs  % can be used to inform pair history

% N trial length items
completed_trials   = behS.n_completed_trials;
behS.hit_history   = beh_dat.saved.PWM_hit_history;             % NaN if violated, 1 if hit, 0 if miss
behS.pair_history  = beh_dat.saved.PWM_pair_history';           % transpose
behS.prev_side     = beh_dat.saved.SideSection_previous_sides'; % transpose
behS.delay         = beh_dat.saved_history.SideSection_Del_time(1:completed_trials);
behS.correct_side  = beh_dat.saved_history.SideSection_ThisTrial(1:completed_trials);
behS.aud1_sigma    = beh_dat.saved_history.PWMSection_AUD1_sigma(1:completed_trials);
behS.aud2_sigma    = beh_dat.saved_history.PWMSection_AUD2_sigma(1:completed_trials);
 
% make a struct w/ rows = trials and columns = pokes, waves, states for
% trial
 for ntrial = 1:completed_trials
     behS.parsed_events(ntrial) = beh_dat.saved_history.ProtocolsSection_parsed_events{ntrial,1}
 end
 
sprintf('behavior info for %s has been extracted', sess_name)
save(save_name,'behS');


