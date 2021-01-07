function find_wireless_sess_E_J(expmtr, rat_name, varargin)
%% ---------------------
% written by Emily Jane Dennis 20200831 - during covid pandemic and
% modified by JRB 20200107 (still in pandemic...)
% purpose is to find the correct bdata file and sessid for ephys
% collections
%
%
% TODO:
% - 
%
% INPUT PARAMETERS:
% - expmtr = experimenter as string
% - ratname = rat name as a string

% OPTIONAL PARAMETERS:
% - dio_path = path to .dio folder(s), assumes cwd otherwise
% - behav_path = path to solo data, assumes `Y:RATTER/SoloData/Data/expmtr`
% - sampling_rate = for ephys, assumes 30000 Hz
% 
% RETURNS:
% - 
% 
% = EXAMPLE CALLS:
% - find_wireless_sess_E_J('Emily, 'W122', 'dio_path', 'W:\jbreda\ephys\W122', 'behav_path', 
% 'Y:\RATTER\SoloData\Data\Emily\')

% ---------------------

%% Assess variable inputs
brody_dir = 'Y:'; % change this depending on how you mount bucket
p = inputParser();
addParameter(p, 'dio_path', '');
addParameter(p, 'behav_path', '');
addParameter(p,'fs',30000);
parse(p,varargin{:});

dio_path = p.Results.dio_path;
behav_path = p.Results.behav_path;
sampling_rate = p.Results.fs;

% if you didn't pass in a dio_path, this function assumes you are running
% it from the rats scratch folder (very specific to JRB kilosort pipeline)
if isempty(dio_path)
    dio_path = fullfile('W:\jbreda\ephys', rat_name);
end

% if you didn't pass in a behavior path, given experimenter & current brody
% lab structure, it will assume the following path
if isempty(behav_path)
    behav_path = fullfile(brody_dir, 'RATTER/SoloData/Data', expmtr);
end
%% Iterate through DIO & match to behavior 
cd(dio_path)
listoffiles = dir([dio_path '/*DIO']); % get list of dio folders

% for each dio in the folder, find session
for i=1:length(listoffiles)
    cd(listoffiles(i).name);
    diofile = dir('*.dat');
    diofilename=diofile.name(1:end-4);
    datafile = readTrodesExtractedDataFile(diofile.name);
    % dio file has both on and off times, we take just the on times and 
    % divide by the sampling rate
    ephys_ttls_sec = datafile.fields(1).data(datafile.fields(2).data==1)/samplingrate;
    totalsec = ephys_ttls_sec(end)-ephys_ttls_sec(1);
    
    % first, pull out the rat (defined by user), date, and length of the recording
    if strcmp(diofile.name(1:4),'data')
        % saved in the format data_sd...
        sessdate = diofile.name(10:17);
    else % saved by ratname_date
        sessdate = [diofile.name(12:15) diofile.name(6:7) diofile.name(9:10)];
    end
    
    % find sessid for the rat and 
    [sessid]=bdata(sprintf('select sessid from sessions where ratname="%s" and sessiondate="%s"',ratname, sessdate));
 
    if isempty(sessid) % if find 0 sessions - return [] and warning
        warning('sessid is empty for rat %s',ratname)
        bdatafile=[];
        minimumvalue=1e10;
    elseif length(sessid)> 1
        % look in bdata file and find timing
        [starttimes] = bdata(sprintf('select starttime from sessions where ratname="%s" and sessiondate="%s"',ratname, sessdate));
        [endtimes] = bdata(sprintf('select endtime from sessions where ratname="%s" and sessiondate="%s"',ratname, sessdate));
        for j=1:length(starttimes)
            sesstiming(j) = seconds(datetime(endtimes{j})-datetime(starttimes{j}))-totalsec;
        end
        [minimumvalue,indexofmin]=min(sesstiming);
        [bdatafile]=bdata(sprintf('select data_file from sessions where sessid="%d"',sessid(indexofmin)));
        
    else     % if find 1 session - great! we're done, pull out sessid
        [bdatafile] = bdata(sprintf('select data_file from sessions where ratname="%s" and sessiondate="%s"',ratname, sessdate));
        minimumvalue=0;
    end
    
    if ~isempty(bdatafile)
        load(sprintf('%s/%s',solo_path,bdatafile{1}))
        bdatasavename = [diofilename '_bdatafile.m'];
        save(bdatasavename,'saved','saved_history')
        save([diofilename 'bdata_ephys_timediff.m'],'minimumvalue')
    end
    
    if python == 1 %let's convert bdata, dio to python
        pyfoldername = [diofilename '_python'];
        mkdir (fullfile(pwd, pyfoldername))
        cd(pyfoldername)
        writeNPY(ephys_ttls_sec,[diofilename '_ephys_ttls'])
        
        for k=1:length(saved_history.ProtocolsSection_parsed_events)
            scopeon(k) = saved_history.ProtocolsSection_parsed_events{k,1}.waves.TrigScope(1,1);
        end
        scopeon=scopeon(2:end);
        % note this could be longer or shorter than the actual ephys ttls
        % depending on when sess was started, if a session was restarted,if
        % batteries failed, or  if testing of the headstage occurred
        writeNPY(scopeon,[diofilename 'bdata_ttls'])
        bdata_to_npy(solo_path, bdatafile{1})
    end
    
    % return to starting folder
    cd(dio_path)

end
        
end