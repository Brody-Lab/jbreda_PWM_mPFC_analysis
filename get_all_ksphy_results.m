function get_all_ksphy_results(sess_path)
 %% ---------------------
% Written by Jess Breda on 2020-01-19 to allow for get_ks_phyresults_J.m
% to be run on multiple sessions. 
%
% INPUT PARAMETERS:
% - sessions_path = string of directory path that has session(s) to be
% found/aligned in bdata as folder name
%
% RETURNS:
% - saves a SpkS structure (see get_ksphy_results for more info)for each 
% session along with cluster notes

% EXAMPLE CALL:
% find_all_wireless_sess('Y:\jbreda\ephys\post_sort_analysis\sorted_pre_bdata')
%
% ---------------------

% get directory infomation
sess_path_info = dir(sess_path);

% grab all the names (aka session names)
all_session_names = {sess_path_info.name};

for i = 3:length(all_session_names) % starting at 3 to skip ". and .."
    % grab the session name & pass it into the find_wireless_sess fx
    session_name = all_session_names{i}
    get_ksphy_results_J(session_name)
end
end
