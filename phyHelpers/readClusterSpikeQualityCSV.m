function [cids, csq] = readSpikeQualityCSV(filename)
% cids is length nClusters, the cluster ID numbers
% csq is length nClusters with spike quality:
% This is specific to JRB curation
% - 1 = m = mua
% - 2 = s = good
fid = fopen(filename);
C = textscan(fid, '%s%s'); 
fclose(fid);

cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, cids);
cids = [cids{~ise}];

isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
isMUA = cellfun(@(x)strcmp(x,'m'),C{2}(2:end));
isGood = cellfun(@(x)strcmp(x,'s'),C{2}(2:end));
csq = zeros(size(cids));

csq(isMUA) = 1;
csq(isGood) = 2;
csq(isUns) = 3;