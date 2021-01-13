function spkS = get_ksphy_results_T(sess, varargin)
% function get_ksphy_waveforms(sess_name)
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
% sync timestamps of timing variable

%% PARSE INPUTS
p = inputParser();
% parent directory w/ folders containing sorted sessions that each have sub
% bundle binary + kilosort output
addParameter(p, 'sess_dir', ''); 
addParameter(p, 'overwrite', false); % whether to reload or overwrite results
addParameter(p, 'rematch', false); % whether to reload or overwrite session match
addParameter(p, 'curator_name', 'Jess'); % name of person curating data
parse(p,varargin{:});

sorted_dir    = p.Results.sess_dir;
curator     = p.Results.curator_name;
overwrite   = p.Results.overwrite;
rematch     = p.Results.rematch;

if ~ispc
   error('the paths for this fx are all specified for PC, change them and then comment this out')
end

%% LOCATE PATHS & FIND BEHAVIOR SESS
% run find wireless session to get directory infomation, these are all
% defuallt inputs, but incase you want to change them i've written it out
sess_match = find_wireless_sess_J(sess, 'overwrite', '0', 'rat_name', 'W122', ...
    'expmtr', 'Emily', 'behav_dir', 'Y:\RATTER\SoloData\Data\W122', ...
    'mdas_dir', 'W:\jbreda\ephys\W122')

% grab mda path from output
save_path = sess_match.res.save_path
[mda_dir, ~] = fileparts(save_path)

% find kilosort info based on input or session name
if isempty(sorted_dir)
    sorted_dir_JB = 'Y:\jbreda\ephys\post_sort_analysis\find_wireless_sess_test\get_kphys_results_test'
    sorted_sess_dir = fullfile(sorted_dir_JB, sess)
else
    sorted_sess_dir = fullfile(sorted_dir, sess)
end
%% BUNDLE & CHANNEL INFO
% I save my files as _firstbundle, _secondbundle, Tyler does _bundle1
% _bundle2 so I need to iterate differntly than his code

% bundles
bndl_names = {'first' 'fourth' 'second' 'third'}
bndl_dirfun = @(n_bndl) fullfile(sorted_sss_dir, sprintf('%s_%sbundle',sess_name,n_bndl));

% channels
mda_filefun = @(n_chan) fullfile(mda_dir, sprintf('%s.nt%i.mda',sess_name,n_chan));

% where to save info & overwrite checl
save_name   = fullfile(sorted_sess_dir,'ksphy_clusters.mat');
if exist(save_name,'file') && ~overwrite
    load(save_name,'spkS');
    return
end
%%
uv_per_bit  = 1;
warning('uv per bit conversion ratio unknown')
nchperb     = 32;
nbundles    = 4; % I might not always have nbundles should make this flexible
wave_x      = -6:25;
nwaves      = 10000;
ch2tt       = @(ch, bb) ceil(ch / 4) + (bb-1)*nchperb/4;
% Loop over bundles to get cluster information and figure out which
% tetrodes we need to load to get cluster waveforms
%S(nbundles) = struct();

notes_path  = fullfile(sorted_sess_dir,'cluster_notes.txt');
notes_fid   = fopen(notes_path,'w+');
ratname     = sess_match.ratname;
sessiondate = datestr(['20' sess_match.date_str([1 2]) '-' sess_match.date_str([3 4]) ...
    '-' sess_match.date_str([5 6])],29);



fprintf(notes_fid,'%s\n%s\n%s\n\n',sessiondate, ratname, curator);


for bb = 1:nbundles;
    % load up spike times and cluster ids from phy
    bundle_dir  = bndl_dirfun(bb);
    cinf_path  = fullfile(bundle_dir,'cluster_info.tsv');
    if ~exist(cinf_path)
        in = input(sprintf('couldn''t find anything for bundle %i. want to continue?', bb),'s');
        if lower(in) == 'y'
            continue;
        else
            error();
        end
    end
    sp          = loadKSdir(bundle_dir);
    sp.mua      = sp.cgs == 1;
    sp.single   = sp.cgs == 2;

    % Find which tetrode each cluster is on using cluster info file
    if ~exist(cinf_path,'file')
        prompt = sprintf(['could not find cluster info file for bundle %i. '...
            'Do you want to continue? (y/n)'],bb);
        in = input(prompt, 's');
        if lower(in) == 'y'
            continue
        else
            return
        end
    end
    fid = fopen(cinf_path);
    C   = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s');
    assert(~isempty(strfind(C{1}{1}, 'id'))); % these ids match the phy gui
    
    cinfo_id = cellfun(@str2num,C{1}(2:end));
    assert(strcmp(C{6}(1), 'ch')); % channel w/ strongest template
    
    cinf_ch = cellfun(@str2num, C{6}(2:end));
    ncids   = length(sp.cids);
    sp.ch1  = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    sp.tt1  = nan(size(sp.cids)); % 1 indexed tetrode id
    sp.cid1 = nan(size(sp.cgs)); % 1 indexed cluster id. numbers span bundles so we don't have non-unique cluster ids
    sp.nspk = nan(size(sp.cgs)); % how many spikes are in each cluster

    % loop over good/mua clusters
    for cc = 1:ncids
        clu_ix      = sp.cids(cc); % get phy cluster id
        info_ix     = cinfo_id == clu_ix; % find index for this cluster in info file
        sp.ch1(cc)  = cinf_ch(info_ix) + 1; % best channel for cluster (indexed from 1)
        sp.tt1(cc)  = ch2tt(sp.ch1(cc),bb); % convert best channel num to best tetrode num
        sp.nspk(cc) = sum(sp.clu == clu_ix);   
    end
    S(bb) = sp;
end
clear sp

% create a filter for the waveforms
assert(sum(diff([S.sample_rate]))==0) % check that all the files have same sampling rate
fs          = S(1).sample_rate;
[n,f0,a0,w] = firpmord([0 1000 6000 6500]/(fs/2), [0 1 0.1], [0.01 0.06 0.01]);
spk_filt    = firpm(n,f0,a0,w,{20});


% Compute each cluster's mean waveform by loading relevant tetrode,
% filtering the signal and pulling relevant indices
nactivetts  = length(unique([S.tt1]));
spkS(nactivetts) = struct();

tt_ix = 0;
for bb = 1:length(S)
    if isempty(S(bb))
        warning(fprintf('skipping bundle %i', bb));
        continue; 
    end
    active_tts = unique(S(bb).tt1);
    fs = S(bb).sample_rate;
    % loop over trodes in this bundle with clusters
    for tt = 1:length(active_tts)
        
        % Figure out which clusters are on this tetrode
        tt_ix        = tt_ix + 1;
        this_tt     = active_tts(tt);
        this_mda    = mda_filefun(this_tt);
        this_tt_ind = S(bb).tt1 == this_tt;
        active_clu  = S(bb).cids(this_tt_ind);
        is_mua      = S(bb).mua(this_tt_ind);
        is_single   = S(bb).single(this_tt_ind);
        
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
        tt_nspk  = sum(S(bb).nspk(this_tt_ind));
        ev_st    = nan(tt_nspk,1);
        ev_ind   = nan(tt_nspk,1);
        tt_clu   = nan(tt_nspk,1);
        
        nchpertt = 4;
        n_clu_on_tt = length(active_clu);
        event_waves = nan(n_clu_on_tt, nchpertt, length(wave_x), nwaves);
        spk_ix_keep = nan(n_clu_on_tt,nwaves); % indices used in mean waveform

        end_ind     = 0; % keep track of last entry into tetrode
        for cc = 1:n_clu_on_tt 
            if is_mua(cc)
                clus_label = sprintf('%i multi\n',cc);
            elseif is_single(cc)
                clus_label = sprintf('%i single\n',cc);
            end
            fprintf(notes_fid, clus_label);
            
            this_cid    = active_clu(cc);
            this_st     = S(bb).st(S(bb).clu == this_cid);
            this_spk_ix = round(this_st * fs);
            this_nspk   = length(this_st);
            start_ind   = end_ind + 1;
            end_ind     = end_ind + this_nspk;
            ind         = start_ind:end_ind;
            ev_st(ind)  = this_st;
            ev_ind(ind) = this_spk_ix;
            tt_clu(ind) = ones(size(this_st)) * cc;
            keep        = 1:min([nwaves this_nspk]);
            rp_spk_ix   = this_spk_ix(randperm(this_nspk));
            spk_ix_keep(cc,keep) = sort(rp_spk_ix(keep));

            for ss = 1:length(keep)
                tmpWf = dat_filt(:,spk_ix_keep(cc,ss)+wave_x);
                event_waves(cc,:,:,ss) = tmpWf;
            end
            
            
        end

        spkS(tt_ix).ratname      = ratname;
        spkS(tt_ix).mua          = is_mua;
        spkS(tt_ix).single       = is_single;
        spkS(tt_ix).recpath      = this_mda;
        spkS(tt_ix).trodenum     = this_tt;
        spkS(tt_ix).event_ind    = ev_ind;
        spkS(tt_ix).event_ts     = ev_st;
        spkS(tt_ix).event_clus   = tt_clu;
        spkS(tt_ix).phy_cids     = active_clu;
        spkS(tt_ix).fs           = fs;
        spkS(tt_ix).event_wave   = -event_waves;
        spkS(tt_ix).wave_x       = wave_x;
        spkS(tt_ix).wave_t_s     = wave_x/fs;
        spkS(tt_ix).waves_mn     = -nanmean(event_waves,4);
        spkS(tt_ix).waves_std    = -nanstd(event_waves,[],4);
        spkS(tt_ix).waves_clus   = 1:n_clu_on_tt;
        spkS(tt_ix).waves_ind    = spk_ix_keep;
        spkS(tt_ix).sess_match   = sess_match;
        spkS(tt_ix).sync_fit_m   = sess_match.spk2fsm_rt(1);
        spkS(tt_ix).sync_fit_b   = sess_match.spk2fsm_rt(2);
        spkS(tt_ix).event_ts_fsm = sess_match.spk2fsm_fn(ev_st);
        spkS(tt_ix).clusnotespath = notes_path;

    end
end

fclose(notes_fid);

save(save_name,'spkS');

