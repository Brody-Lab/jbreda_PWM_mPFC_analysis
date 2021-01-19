function spkS = get_ksphy_results_J(sess, varargin)
% function get_ksphy_waveforms(sess_name)
%
% Written by Tyler Boyd-Meridith and adjusted by Jess Breda on 2021-01-15
% for analysis of wireless tetrode data
% -----Tyler's notes-----
% For each bundle
    % 1. use sp = loadKSdir to get the spike times and cluster ids of all mua and good clusters
        % convert sp.st into an index into the mda file (might be some
        % annoying offsets here?)
    % 2. use the cluster_info.tsv file to figure out which tetrode each cluster is
    % on
    % 3. For each tetrode
        % 1. load the relevant mda file (make sure its in uv)
        % 2. grab timepoints around each of N spikes
        % 3. average each clusters spikes together
        % 4. Append the spike times, cluster labels (relabeled to be
        % continuous + separate struct of phy labels), and average
        % waveforms and sem of waveforms - i could also try to save all
        % waveforms but i think these files will get too big - can make it
        % an optional argument

% get the syncing parameters
% sync timestamps of timing variable %maybe add figure path to find_wireles
% -------JRB--------
% INPUT PARAMETERS:
% - sess = name of sorted session you'd like to align & get info for
%
% OPTIONAL PARAMETERS:
% - sorted_dir =  parent directory w/ folders containing sorted sessions 
% each have sub folder sof reach bundle with kilsort output
% - overwrite = T/F whether to overwrite results of this fx if already run
% - curator_name = name of person curating the data as a string
% !!NOTE!! see call for find_wireless_sess_J to adjust input to that fx as
% needed
%
% RETURNS:
% 
% 
% EXAMPLE CALL:
% 
%
% TODO:
% - line 182 onward
%
% assumes you are running on a PC!


%% PARSE INPUTS
sess = 'data_sdc_20190902_145404_fromSD' % for debugging comment out later

p = inputParser();
addParameter(p, 'sorted_dir', ''); 
addParameter(p, 'overwrite', true); 
addParameter(p, 'curator_name', 'Jess'); 
parse(p,varargin{:});

sorted_dir  = p.Results.sorted_dir;
curator     = p.Results.curator_name;
overwrite   = p.Results.overwrite;

if ~ispc
   error('the paths for this fx are all specified for PC, change them and then comment this out')
end

%% LOCATE PATHS & FIND BEHAVIOR SESS
% run find wireless session to get directory infomation, these are all
% defuallt inputs, but incase you want to change them i've written it out
sess_match = find_wireless_sess_J(sess, 'overwrite', 0, 'rat_name', 'W122', ...
    'expmtr', 'Emily', 'behav_dir', 'Y:\RATTER\SoloData\Data\Emily', ...
    'mdas_dir', 'W:\jbreda\ephys\W122')

% grab mda path from output
save_path = sess_match.save_path
[mda_dir, ~] = fileparts(save_path)

% find kilosort info based on input or session name
if isempty(sorted_dir)
    sorted_dir_JB = 'Y:\jbreda\ephys\post_sort_analysis\get_kphys_results_test';
    sorted_sess_dir = fullfile(sorted_dir_JB, sess);
else
    sorted_sess_dir = fullfile(sorted_dir, sess);
end
%% BUNDLE, CHANNEL & SAVING INFO

% global vars
uv_per_bit  = 1;
warning('uv per bit conversion ratio unknown')
nchperb     = 32;
nbundles    = 4; 

% make bundle dir (stored as sorted_sess_dir\sess_name_bundle1_forkilossort)
bndl_dirfun = @(n_bndl) fullfile(sorted_sess_dir, sprintf('%s_bundle%i_forkilosort',sess,n_bndl));

% make mda channel dir (sorted as mda_dir/session_name_nt1.mda)
mda_filefun = @(n_chan) fullfile(mda_dir, sprintf('%s.nt%i.mda',sess,n_chan));

% go from phy channels to tetrode channles (bc running separate bundles)
ch2tt       = @(ch, n_bndl) ceil(ch / 4) + (n_bndl-1)*nchperb/4;

% bdata save out info
wave_x      = -6:25; % n time pts rel. to each spike that get included in the bdata cells and spktimes database
nwaves      = 10000; % how many waveforms to store/use to compute the average waveform

% where to save info & overwrite check
save_name  = fullfile(sorted_sess_dir,'ksphy_clusters.mat');
save_name2 = fullfile(sorted_sess_dir,'ksphy_clusters_JB.mat');
if exist(save_name,'file') && ~overwrite
    load(save_name,'spkS');
    return
end

% where to save cluster_notes file & adding base info to it
notes_path  = fullfile(sorted_sess_dir,'cluster_notes.txt');
notes_fid   = fopen(notes_path,'w+');
ratname     = sess_match.ratname;
sessiondate = datestr(['20' sess_match.date_str([1 2]) '-' sess_match.date_str([3 4]) ...
    '-' sess_match.date_str([5 6])],29);

fprintf(notes_fid,'%s\n%s\n%s\n\n',sessiondate, ratname, curator);

%% INTERACT W/ PHY & FIND CHANNELS WITH CELLS

% Loop over bundles to get cluster information and figure out which
% tetrodes we need to load to get cluster waveforms
for n_bndl = 1:nbundles;
    
    % load up spike times and cluster ids from phy
    bundle_dir  = bndl_dirfun(n_bndl);
    cinf_path  = fullfile(bundle_dir,'cluster_info.tsv');
    
    
    % Phy helper will load spike quality/cluster group info per JRB edits
    % the function currently exludes any templates marked as Noise
    sp = loadKSdir(bundle_dir); 
    
    % if there is an spike quality file, use it to assign MUA/Single
    % NOTE: specific to how JRB sorts!
    if isfield(sp, 'csq') 
        sp.mua      = sp.csq == 1; 
        sp.single   = sp.csq == 2;
        sp.sort_metric = 'spike_quality'
    else
        sp.mua      = sp.cgs == 1;
        sp.single  = sp.cgs == 2;
        sp.csq = [] % need to do this to keep structure size stable
        sp.sort_metric = 'cluster_group'
    end
    
    if ~exist(cinf_path)
        sprintf('couldn''t find any cluster_info for bundle %i, will continue.', n_bndl)   
    end

    % access cluster_info.tsv & extract for single/muas
    fid = fopen(cinf_path);
    C   = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s'); % JRB has extra string here bc dded a 'sq' column
    
    assert(~isempty(strfind(C{1}{1}, 'id'))); % these ids match the phy gui, make sure they are there
    cinfo_id = cellfun(@str2num,C{1}(2:end)); % get ids
    
    assert(strcmp(C{6}(1), 'ch'));            % channel w/ strongest template
    cinf_ch = cellfun(@str2num, C{6}(2:end)); % get strongest templates 
    
    % initialize space
    ncids   = length(sp.cids);    % number of good/muas present
    sp.ch1  = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    sp.tt1  = nan(size(sp.cids)); % 1 indexed tetrode id
    sp.cid1 = nan(size(sp.cgs));  % 1 indexed cluster id. numbers span bundles so we don't have non-unique cluster ids
    sp.nspk = nan(size(sp.cgs));  % how many spikes are in each cluster

    % loop over good/mua clusters
    for cc = 1:ncids  
        clu_ix      = sp.cids(cc);              % get phy cluster id
        info_ix     = cinfo_id == clu_ix;       % find index for this cluster in info file
        sp.ch1(cc)  = cinf_ch(info_ix) + 1;     % best channel for cluster (indexed from 1)
        sp.tt1(cc)  = ch2tt(sp.ch1(cc),n_bndl); % convert best channel num to best tetrode num
        sp.nspk(cc) = sum(sp.clu == clu_ix);    % get number of spikes for that cluster  
    end
    S(n_bndl) = sp;                             % save out
end
clear sp

%% GET SPIKE TIMEs & IDXS
% create a filter for the waveforms
assert(sum(diff([S.sample_rate]))==0) % check that all the files have same sampling rate
fs          = S(1).sample_rate;
              %firpmord(cut off freq, dsired amplitudes, maximum allowed deviation           
[n,f0,a0,w] = firpmord([0 1000 6000 6500]/(fs/2), [0 1 0.1], [0.01 0.06 0.01]);
              %firmp(filter order, normalized freq pts, desired amplitude, weights, ?) 
spk_filt    = firpm(n,f0,a0,w,{20});

% Compute each cluster's mean waveform by loading relevant tetrode,
% filtering the signal and pulling relevant indices

%initialize
nactivetts  = length(unique([S.tt1]));
spkS(nactivetts) = struct();
tt_ix = 0;

% JB non SpkS storage
n_clu = length([S.tt1]);
PWMspkS(n_clu) = struct();
clu_ix = 0;

for n_bndl = 1:length(S)
    if isempty(S(n_bndl))
        warning(fprintf('skipping bundle %i', n_bndl));
        continue; 
    end
    active_tts = unique(S(n_bndl).tt1);
    fs = S(n_bndl).sample_rate;
    
    % loop over channels in this bundle with clusters
    for n_chan = 1:length(active_tts)
        
        % Figure out which clusters are on this tetrode
        tt_ix        = tt_ix + 1;                   
        this_tt     = active_tts(n_chan);            % find active channel    
        this_mda    = mda_filefun(this_tt);          % get mda file directory
        this_tt_ind = S(n_bndl).tt1 == this_tt;      % create array of same size as number of clusters on chan
        active_clu  = S(n_bndl).cids(this_tt_ind);
        is_mua      = S(n_bndl).mua(this_tt_ind);    % determine cluster type
        is_single   = S(n_bndl).single(this_tt_ind);
        
        fprintf(notes_fid, "TT%i - \n",this_tt); 
        
        % Load this tetrode and filter it
        fprintf('loading mda file for tetrode %i...',this_tt);
        tic;
        dat    = readmda(this_mda);
        toc;
        fprintf('filtering');
        tic;
        dat_filt    = filtfilt(spk_filt,1,uv_per_bit*dat')';
        toc
        
        % intialize variables of interest for this tetrode
        tt_nspk  = sum(S(n_bndl).nspk(this_tt_ind));
        ev_st    = nan(tt_nspk,1);   
        ev_ind   = nan(tt_nspk,1);
        tt_clu   = nan(tt_nspk,1);
        
        nchpertt = 4;
        n_clu_on_tt = length(active_clu);
        event_waves = nan(n_clu_on_tt, nchpertt, length(wave_x), nwaves);
        spk_ix_keep = nan(n_clu_on_tt,nwaves); % indices used in mean waveform

        end_ind     = 0; % keep track of last entry into tetrode to add spike info to strcut
     
        for cc = 1:n_clu_on_tt 
            clu_ix = clu_ix + 1
            if is_mua(cc)
                clus_label = sprintf('%i multi\n',cc);
            elseif is_single(cc)
                clus_label = sprintf('%i single\n',cc);
            end
            fprintf(notes_fid, clus_label);
            
            this_cid    = active_clu(cc);                          % cluster id
            this_st     = S(n_bndl).st(S(n_bndl).clu == this_cid); % cluster spike times in s
            this_spk_ix = round(this_st * fs);                     % idx value for each spike
            this_nspk   = length(this_st);                         % n spikes
            start_ind   = end_ind + 1;                             % create idx for varying size of spikes per clus/tetrode
            end_ind     = end_ind + this_nspk;
            ind         = start_ind:end_ind;                        
            ev_st(ind)  = this_st;                                 % store spike times                          
            ev_ind(ind) = this_spk_ix;                             % store spike idx
            tt_clu(ind) = ones(size(this_st)) * cc;                 
            keep        = 1:min([nwaves this_nspk]);               % time window size to keep for avg waveform
            rp_spk_ix   = this_spk_ix(randperm(this_nspk));        % picking 10000 rand spike idx to keep
            spk_ix_keep(cc,keep) = sort(rp_spk_ix(keep));

            for ss = 1:length(keep)                                % for each rand spk idx, grab a window of size wave & store                          
                tmpWf = dat_filt(:,spk_ix_keep(cc,ss)+wave_x);
                event_waves(cc,:,:,ss) = tmpWf;
            end
            
            PWMspkS(clu_ix).ratname = ratname;         % rat name  
            PWMspkS(clu_ix).date = sess_match.date_str % ephys session date
            PWMspkS(clu_ix).sessid = sess_match.sessid % behav session id
            PWMspkS(clu_ix).recpath = this_mda;
            PWMspkS(clu_ix).mua          = is_mua;     % is multi?
            PWMspkS(cly_ix).single       = is_single;
            PWMspkS(clu_ix).trodenum = this_tt;
            PWMspkS(clu_ix).phy_cids  = active_clu;
            PWMspkS(clu_ix).fs           = fs; 
            PWMspkS(clu_ix).event_ind = this_spk_ix;
            PWMspkS(clu_ix).ecent_ts = this_st;
            PWMspkS(clu_ix).sync_fit_m   = sess_match.spk2fsm_rt(1);     % alignment slope
            PWMspkS(clu_ix).sync_fit_b   = sess_match.spk2fsm_rt(2);     % alignment y int
            PWMspkS(clu_ix).event_ts_fsm = sess_match.spk2fsm_fn(this_st); % spike times in fsm time in seconds
            PWMspkS(clu_ix).clusnotespath = notes_path;  
            
            
        end

        spkS(tt_ix).ratname      = ratname;    % rat name                
        spkS(tt_ix).mua          = is_mua;     % is multi?
        spkS(tt_ix).single       = is_single;  % is single?
        spkS(tt_ix).recpath      = this_mda;   % path to mda
        spkS(tt_ix).trodenum     = this_tt;    % tetrode number (1-32)
        spkS(tt_ix).event_ind    = ev_ind;     % spike time indices 
        spkS(tt_ix).event_ts     = ev_st;      % spike time in seconds
        spkS(tt_ix).event_clus   = tt_clu;     % idk
        spkS(tt_ix).phy_cids     = active_clu; % phy cluster ids
        spkS(tt_ix).fs           = fs;         % sampling rate
        spkS(tt_ix).event_wave   = -event_waves;                 % idk
        spkS(tt_ix).wave_x       = wave_x;                       % n time point relative to spike included
        spkS(tt_ix).wave_t_s     = wave_x/fs;                    % wave_x in seconds
        spkS(tt_ix).waves_mn     = -nanmean(event_waves,4);      % mean waveform in uv
        spkS(tt_ix).waves_std    = -nanstd(event_waves,[],4);    % std of waveform
        spkS(tt_ix).waves_clus   = 1:n_clu_on_tt;                % matches waves_mn to event_clus and cluster notes
        spkS(tt_ix).waves_ind    = spk_ix_keep;                  % indices of waves used to compute mean 
        spkS(tt_ix).sess_match   = sess_match;                   % beahvior alignment info
        spkS(tt_ix).sync_fit_m   = sess_match.spk2fsm_rt(1);     % alignment slope
        spkS(tt_ix).sync_fit_b   = sess_match.spk2fsm_rt(2);     % alignment y int
        spkS(tt_ix).event_ts_fsm = sess_match.spk2fsm_fn(ev_st); % spike times in fsm time in seconds
        spkS(tt_ix).clusnotespath = notes_path;                  % path to notes

    end
end

fclose(notes_fid);

save(save_name,'spkS');
save(save_name2,'PWMspkS');

