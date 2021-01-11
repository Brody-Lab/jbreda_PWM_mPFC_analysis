function res = find_wireless_sess_J(sess, varargin)
% %% ---------------------
% orifinally written by Tyler Boyd-Meridth and adapted by Jess Breda
% 20210111. The purpose of this function is to find the correct bdata file
% & session id for ephys collection & save out information for ttl alignmnet
% between finite state machine (fsm) and trodes ephys
%
% INPUT PARAMETERS:
% - sess = name of session you'd like to find & align as a string.
% - rat_name = name of rat in session as a string 
%
% OPTIONAL PARAMETERS:
% - overwrite = if this session has already been run, 0 to skip, 1 to
% rerun
% - expmtr = experimenter who ran the physiology as a string 
% - behav_dir = directory where behavior data is stored e.g. 'Y:\RATTER\SoloData\Data\Emily\'
% - mdas_dir = directory where .mda/.dio folders are for this (and other)
% sessions
% - fs = sampling rate of electrophsyiology recording 
% 
% RETURNS:
% - if match is found for a session, saves the file ttl_match.mat in the 
% session.mda folder with the session & rat name passed into the function 
% as well as a res struct with
% paths, ttl conversions & additonal info needed for alignment.
% 
% EXAMPLE CALL:
% find_wireless_sess_J('data_sdb_20190724_193007_fromSD', 'W122', 'expmtr', 'Emily', 'behav_dir', 
% 'Y:\RATTER\SoloData\Data\Emily\')
%
% TODO:
% - rat_name iteration fix
% to iterate over .mda folders or not?
%
% assumes you are running on a PC!
%
% ---------------------
%% PARSE INPUTS
p = inputParser();
addParameter(p,'overwrite',1); % defualt is skip if alrady run
addParameter(p,'rat_name', 'W122'); % if multiple rats: {'H191','H176'}
addParameter(p,'expmtr','Emily');
addParameter(p,'behav_dir',''); % will assign below
addParameter(p,'mdas_dir',''); % will assign below
addParameter(p,'fs',30000)
parse(p,varargin{:});

brody_dir = 'Y:'; % change this depending on how you mount bucket
overwrite   = p.Results.overwrite;
rat_name     = p.Results.rat_name;
expmtr      = p.Results.expmtr;
behav_dir   = p.Results.behav_dir;
mdas_dir     = p.Results.mdas_dir;
fs          = p.Results.fs;

%% LOCATE PATHS 
% --- behavior data --- 
% if there is no behavior directory passed in, use the experimenters
% behavior directory 
if isempty(behav_dir)
    behav_dir = fullfile(brody_dir, 'RATTER/SoloData/Data', expmtr);
end

% --- ephys ---
% if there is no mda directory, use my folder where I store ephys by rat 
% name, look in rat_names folder, and find the .mda folder for the
% session (same thing for if mdas_dir was passed in)
if  isempty(mdas_dir)
    ephys_dir    = fullfile('W:\jbreda\ephys\', rat_name);
    mda_dir     = fullfile(ephys_dir, sprintf('%s.mda', sess));
else 
    mda_dir = fullfile(mdas_dir, sprintf('%s.mda', sess));
end

% make sure the mda_dir is a folder & create a .mat file to save session & 
% alignment into later. in get_ksphy_results, you want this to be in .mda
% folder
assert(exist(mda_dir,'dir')==7) % 7 = name is a folder
save_path   = fullfile(mda_dir, 'ttl_match.mat');

% if this file has already been run, check overwrite status
if exist(save_path) & ~overwrite
    dio = load(save_path,'res','sess');
    if strcmp(dio.sess, sess)
        res = dio.res;
        if isfield(res, 'spk2fsm_rt')
            return;
        end
    else
        clear dio
    end
end

% use .mda folder name to create dio file name & make sure it is there
% current name structur is: sess.DIO/data_sdb_20190724_193007_fromSD.dio_RFINPUT.dat)
dio_file    = [mda_dir sess '.dio_RFINPUT.dat'];
if ~exist(dio_file)
    % repace the .mda with .DIO to change into the correct directory & grab
    % the correct fname 
    dio_file = fullfile([strrep(mda_dir,'.mda','.DIO')], [sess '.dio_RFINPUT.dat']);
end

fi=dir(dio_file);
if(isempty(fi))
    error(sprintf('dio missing %s',dio_file))
end

%% GRAB DATE
% sometimes sessions are named like 1) data_sdc_20190730_143641_fromSD or
% 2) W122_09_21_2019_1_fromSD and this code will extract the meaningful
% info from each title & reformat into a usable data for bdata search

if strcmp(sess(1:4),'data')
    is_yyyy = regexp(sess,'20\d\d\d\d\d\d_\d\d');
    if ~isempty(is_yyyy)
        ind = is_yyyy + 2;
    else
        ind = regexp(sess,'\d\d\d\d\d\d_\d\d');
    end
    date_str = sess(ind + (0:5));
    fprintf('name in data_YYMMDD format; session data is: %s',date_str)
    
else
    date_str = sess([14 15 6 7 9 10]);
    fprintf('named in rat_MM_DD_YYYY format; session date is: %s',date_str)
end

% get date in SoloData format
sessiondate = ['20' date_str([1 2]) '-' date_str([3 4]) ...
    '-' date_str([5 6])];

%% STOPPED HERE
% load ttls from wireless
dio     = readTrodesExtractedDataFile(dio_file);
c       = dio.fields(1).data;
b       = dio.fields(2).data;
ttls    = double(c(b==1))/fs; % this is ttls in seconds

inc     = 0;

ndays = 60;

for ss = 0:(ndays-1);
    this_datestr = datestr(datenum(sessiondate) - ss,'yymmdd');
    
    for rr = 1:length(rat_name)
        fntemp      = fullfile(behav_dir, rat_name{rr}, ['data_*' this_datestr '*.mat']);
        ratfiles    = dir(fntemp);
        if isempty(ratfiles)
            fprintf(['couldn''t find a match for ' strrep(fntemp,'\','\\')])
            continue
        end
        for ff=1:length(ratfiles)
            inc             = inc +1;
            rats{inc}       = rat_name{rr};
            this_behav_name = ratfiles(ff).name;
            this_behav_path = fullfile(ratfiles(ff).folder, this_behav_name);
            filename{inc}   = this_behav_name;
            behav_path{inc} = this_behav_path;
            warning('off')
            load(this_behav_path, 'saved_history');
            warning('on')
            parsed_events   = saved_history.ProtocolsSection_parsed_events;
            % load behavior TTLS for alignment with ephys
            ttls_fsm = nan(1,length(parsed_events));
            for i=1:length(parsed_events)
                if(isfield(parsed_events{i}.states,'sending_trialnum') && ...
                        length(parsed_events{i}.states.sending_trialnum)>0 )
                    ttls_fsm(i) = parsed_events{i}.states.sending_trialnum(1);
                else
                    ttls_fsm(i) = NaN;
                end
            end
            ttls_fsm=ttls_fsm';
            
            % get interpulse interval of ttls
            trode_gaps  = diff(ttls);
            fsm_gaps    = diff(ttls_fsm);
            
            % look only at ttls delivered 5-10 seconds apart
            trode_gaps_valid    = find(trode_gaps > 5 & trode_gaps < 10);
            % look for correspondance between the phys and behav sess ITIs
            ttls_fsmall=[];
            ttl_trodesall=[];
            for i=1:length(trode_gaps_valid)
                fsm_trodes_diffs=find(abs(fsm_gaps - trode_gaps(trode_gaps_valid(i))) < 0.05);
                for k=1:length(fsm_trodes_diffs)
                    indb=fsm_trodes_diffs(k);
                    if (indb+5)>length(fsm_gaps)
                        continue
                    end
                    if (trode_gaps_valid(i)+5)>length(trode_gaps)
                        continue
                    end
                    vec1    = trode_gaps(trode_gaps_valid(i)+(0:5));
                    vec2    = fsm_gaps(indb+(0:5));
                    
                    %if one of the distances is too long (timeout), abort
                    if(max(abs(diff(vec1)))>30)
                        continue
                    end
                    
                    if(norm(vec1-vec2)>.1)
                        continue
                    else
                        ttls(trode_gaps_valid(i));
                        ttl_trodesall=[ttl_trodesall; ttls(trode_gaps_valid(i)+(0:5))];
                        ttls_fsmall=[ttls_fsmall; ttls_fsm(indb+(0:5))];
                    end
                end
            end
            
            rt{inc}     = polyfit(ttl_trodesall,ttls_fsmall,1);
            
            ttl_trode_fsm = ttl_trodesall*rt{inc}(1)+rt{inc}(2);
            
            if(isempty(ttl_trode_fsm))
                max_residual(inc)   = 100;
                totaldur(inc)       = 0;
            else
                max_residual(inc)   = max(abs(ttl_trode_fsm-ttls_fsmall));
                totaldur(inc)       = (max(ttls_fsmall)-min(ttls_fsmall))/60;
            end
        end
    end
end

% figure it out if any of the behav sessions are acceptable matches
match_ind  = find(max_residual<0.02 & totaldur>5);

if(length(match_ind)<1)
    warning('no match!');
    save_path   = fullfile(mda_dir, 'no_ttl_match.mat');
elseif(length(match_ind)>1)
    disp(valore(match_ind))
    disp(filename(match_ind))
    warning('too many matches!')
    save_path   = fullfile(mda_dir, 'multiple_ttl_matches.mat');
else
    res.goodval     = max_residual(match_ind);
    res.gooddur     = totaldur(match_ind);
    res.spk2fsm_rt  = rt{match_ind};
    res.spk2fsm_fn  = @(spk_ts) spk_ts.*res.spk2fsm_rt(1) + res.spk2fsm_rt(2);
    res.goodpath    = behav_path{match_ind};
    res.behfile     = filename{match_ind};
    res.ratname     = rats{match_ind};
    if bdata('connect') > 0
        sessid          =  bdata(['select sessid from sessions where data_file="'  res.behfile(1:end-4) '"']);
        res.sessid      = sessid;
    end
    res.sessiondate = sessiondate;
end

res.date_str    = date_str;
res.save_path   = save_path;
save(save_path, 'res', 'sess')
