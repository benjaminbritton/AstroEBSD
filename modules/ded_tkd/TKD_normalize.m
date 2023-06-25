function [normalizeddata] = TKD_normalize(fillingdomains,rawdata,numfiles)
%TKD_NORMALIZE Summary of this function goes here
%   Detailed explanation goes here

normalizeddata = 0*rawdata;

% Calculate normalization factors
eventratio = zeros(length(fillingdomains{1}),numfiles-1);
for j=2:numfiles
    reference_data = rawdata(:,1);
    current_data = rawdata(:,j);
    eventratio(:,j-1) = reference_data(fillingdomains{1})./ current_data(fillingdomains{1});
end

eventratio(isinf(eventratio)) = NaN; % account for dead pixels

clear reference_data current_data;

scale_factor = mean(eventratio, 1 ,"omitnan");

for i=1:numfiles
    if i==1
        normalizeddata(:,i) = rawdata(:,i);
    else
        normalizeddata(:,i)= rawdata(:,i) .* scale_factor(i-1);
    end
end
end

