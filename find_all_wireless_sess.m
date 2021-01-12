function find_all_wireless_sess(sess_path)
 %% ---------------------
% Written by Jess Breda on 2020-01-12 to allow for find_wireless_sess_J.m
% to be run on multiple sessions. I am running this on sessions I have
% already sorted and they are stored in a direcoty with session name =
% folder name.
%
% INPUT PARAMETERS:
% - sessions_path = string of directory path that has session(s) to be
% found/aligned in bdata as folder name
%
% RETURNS:
% - - if match is found for a session, saves the file ttl_match.mat in the 
% session.mda folder with alignment information, see find_wireless_sess_J
% for path and more documentation
% 
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
    find_wireless_sess_J(session_name)
end
end
