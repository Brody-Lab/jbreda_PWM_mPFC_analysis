
function get_mask_info(bin_folder, save_folder)

% ---------------------
% adapted from jbreda_kilosort repo kilosort_preprocess_forcluster 
% goal is to get the mask that was applied to .bin files to use for
% analysis (bc i forogot to save it...)
%
% INPUT PARAMETERS:
% - bin_folder folder with .bin files that haven't been processed
% - save_folder: where you want masks sent for saving
% - 
% ---------------------
%% inputs

cd(bin_folder) 

% hard coding some things here 2021-1-26 since they were already run
threshold = 0.5;
window = 10000;
chan=32;
ops.fs     = 32000;    
ops.fshigh = 300;


%% filter & initialization

% make a filter for the data based on ops inputs/defaults 
% butterworth filter with only 3 nodes (otherwise it's unstable for float32)
% Wn (ops.fshigh/ops.fs) is the cutoff frequency, and must be between 0.0 and 1.0. 1.0 is half the
% sample rate
% high means it's a highpass filter
%outputs filter coefficients b1 (numerator) and a1 (denominator)
[b1, a1] = butter(3, ops.fshigh/(ops.fs/2), 'high');
    
% make list of files to process
listofbinaryfiles=dir('*.bin');
%make empty fftosave
fftosave=[];

%determine computer type
if ispc
    delim='\';
else
    delim='/';
end

%% loop over binary files

addpath(save_folder)
for i = 1:length(listofbinaryfiles)

    % first, open the binary file to read
    fname = listofbinaryfiles(i).name;
    sess_name = fname(1:end-4);
    fid=fopen(fname,'r');
    
    full_mask = [];
    cd(save_folder)
    
    
    while 1
    % now, read in a PORTION of the data. Format it as a matrix with chan rows and
    % 1e5 values - we will loop through this for each file until all data is
    % read in
        dataRAW = fread(fid, [chan 1e6], 'int16');
        sizeofdata=size(dataRAW);
        if sizeofdata(2) == 0
            break %breaks the while loop
        end

        % transpose
        dataRAW = dataRAW';
        % divide by 1000 because the filter prefers that
        dataRAW = double(dataRAW)/1000;
               
        % apply the filter
        datr = filtfilt(b1, a1, dataRAW);
        dataFILT = datr; %renaming now so I can plot later without overwriting

        % make a binary 'mask' for the data, looking for big changes to
        % remove
            % first we want absolute values, so we square datr, take the means of each
            % row (channel), and then take the square root
        dataABSVAL = mean(datr.^2, 2).^.5; 
            % create a binary mask where dataABSVAL > threshold
            % then, take moving mean of mask with window size
        mask1MEAN = movmean(double(dataABSVAL>threshold), window);
        
        % binary mask 'signal' = 1, 'noise' = 0
        mask2= mask1MEAN < 0.00001;
        
        % mask data & set noise to 0
        % dataMASK = datr .* mask2;
        
        full_mask = vertcat(full_mask, mask2); % append the masks together
        
        % save
%         dat16 = int16(1000*dataMASK');
%         fwrite(fidw, dat16, 'int16');

    end
    
          
    if contains(sess_name, "first")
        new_name = strcat(sess_name(1:end-12), "_bundle1_mask_info");
          
        
    elseif contains(sess_name, "second")
        new_name = strcat(sess_name(1:end-13), "_bundle2_mask_info");
            
            
    elseif contains(sess_name, "third")
        new_name = strcat(sess_name(1:end-12), "_bundle3_mask_info");
            
            
    elseif contains(sess_name, "fourth")
        new_name = strcat(sess_name(1:end-13), "_bundle4_mask_info");       
            
    end

    writeNPY(full_mask, new_name)
    
    fclose(fid);
    
    cd(bin_folder)

    sprintf('finished file %d of %d files to process',i,length(listofbinaryfiles))
%cd .. %return to directory with other binary files

end
