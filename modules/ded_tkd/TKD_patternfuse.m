function [fusedpattern] = TKD_patternfuse(rawdata,normalizeddata,fillingdomains)
%TKD_PATTERNFUSE Summary of this function goes here
%   Detailed explanation goes here
numfiles=size(rawdata,2);

fusedpattern = zeros(254^2,1);

for i=1:numfiles
    current_data = normalizeddata(:,i);
    fusedpattern(fillingdomains{i}) = current_data(fillingdomains{i});
end

fusedpattern = reshape(fusedpattern, 254,254);
stitched_mask = medfilt2(fusedpattern, [3 3]);
fusedpattern(rawdata(:,1)==0) = stitched_mask(rawdata(:,1)==0);



end

