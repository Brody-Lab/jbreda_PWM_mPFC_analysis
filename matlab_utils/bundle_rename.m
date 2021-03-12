%% FILE RENAME
% using this script to rename files from the format
% 'W122_05_27_2019_1_fromSD_firstbundle_T5_W10000_forkilosort' to
% 'W122_05_27_2019_1_fromSD_bundle1_forkilosort; to better work with tylers
% code

%%
sessions_dir = 'Y:\jbreda\ephys\post_sort_analysis\sorted_pre_bdata'

dir_info = dir(sessions_dir);
sess_names = {dir_info.name};

for i=3:length(sess_names)
    
    % get session folder info
    sess_dir = fullfile(sessions_dir, sess_names{i});
    sorted_dir_info = dir(sess_dir);
    
    % get names of sorted sessions
    sorted_names = {sorted_dir_info.name};
    cd(sess_dir) % for renaming is
    
    % for each sorted file, update the name
    for k=3:length(sorted_names)
    
        old_name = sorted_names{k};
        
        if contains(sorted_names{k}, "first")
            new_name = strcat(sess_names{i}, "_bundle1_forkilosort");
            movefile(old_name, new_name)
        
        elseif contains(sorted_names{k}, "second")
            new_name = strcat(sess_names{i}, "_bundle2_forkilosort");
            movefile(old_name, new_name)
            
        elseif contains(sorted_names{k}, "third")
            new_name = strcat(sess_names{i}, "_bundle3_forkilosort");
            movefile(old_name, new_name)
            
        elseif contains(sorted_names{k}, "fourth")
            new_name = strcat(sess_names{i}, "_bundle4_forkilosort");
            movefile(old_name, new_name)
            
        else
            fprintf('name is: %s and cannot rename', sorted_names{k})
        end
    end
    end
    
  