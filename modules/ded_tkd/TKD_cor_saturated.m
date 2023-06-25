function [fillingdomains, satdomains] = TKD_cor_saturated(rawdata,numfiles,sat_index)
%TKD_COR_SATURATED Summary of this function goes here
%   Detailed explanation goes here
%% Identify saturated domains 
satdomains = cell(1,numfiles);
for i=1:numfiles
    satdomains{i} = find(rawdata(:,i) >= sat_index);
end

fillingdomains = cell(1,numfiles);
for i=1:numfiles
    if i==1
        fillingdomains{i} = find(rawdata(:,i) < sat_index);
    else
        fillingdomains{i} = setdiff(satdomains{i-1}, satdomains{i});
    end
end
end

