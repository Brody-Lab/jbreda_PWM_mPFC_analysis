function find_wireless_sess_E(dio_path,solo_path,ratname,python)
% ---------------------
% written by Emily Jane Dennis 20200831 - during covid pandemic!
% purpose is to find the correct bdata file and sessid for ephys
% collections
%
%
% TODO:
% - 
%
% INPUT PARAMETERS:
% - dio_path = a folder path to .DIO folders
% - solo_path = where the SoloData is stored for that rat e.g. "X:\RATTER\SoloData\Data\Emily\W122"
% - ratename = rat name as a string
% - python = flag == 1 if you want to convert to python
% 
% OPTIONAL PARAMETERS:
% - 
% 
% RETURNS:
% - 
% 
% = EXAMPLE CALLS:
% - find_wireless_sess('path/to/dio/folder(s)', 'path/to/Solo/data;, 'W122', 1)
% ---------------------

%if you have a different sampling rate, change it here
% TODO allow for flexible input with default to 30000
samplingrate = 30000;

% change this for your machine
% TODO make flexible/auto ID
solo_path = "X:\RATTER\SoloData\Data\Emily\W122";

cd(dio_path)
% list of dio folders
listoffiles = dir([dio_path '/*DIO']);

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