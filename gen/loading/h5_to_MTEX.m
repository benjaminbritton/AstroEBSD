function [ebsd] = h5_to_MTEX(Data_InputMap,MTEX_Settings,cslist,phasedata)
%H5_TO_MTEX This will convert a Data_InputMap into an ebsd container
%
% [ebsd] = h5_to_MTEX(Data_InputMap,MTEX_Settings,cslist,phasedata)
% INPUTS:
% MTEX_Settings - structure, 
%                 with optional fields xpref and ypref for axis directions
%
% cslist - crystal symmetry list, 
%          created from e.g. cslist={'notIndexed',loadCIF('Fe-Iron-alpha.cif');}; 
% 
% Data_InputMap - structure with 
%                 phi1, PHI, phi2 - in degrees
%                 XSample, YSample - in micrometers
%                 xpts, ypts = number of points in each dimensions
%
%                 and other properties 
%                 --> these will be read dynamically 
%                     and loaded into ebsd.prop
% 
% If MADPhase exists & phase data is not input, this will be used as the phase input
%
% Outputs
% ebsd - an ebsd container with props and orientations generated


%set preferences for MTEX if needed
if isfield(MTEX_Settings,'xpref')
    setMTEXpref('xAxisDirection',MTEX_Settings.xpref);
end
if isfield(MTEX_Settings,'ypref')
    setMTEXpref('zAxisDirection',MTEX_Settings.ypref);
end

%read pts from the array
xpts=Data_InputMap.xpts;
ypts=Data_InputMap.ypts;

%read the fields
fields_now=fieldnames(Data_InputMap);

%count the fields
numfields=length(fields_now);

%build the prop structure
prop=struct;

%go through each field and determine if it is map data
for n=1:numfields
    %get the data
    data_now=Data_InputMap.(fields_now{n});
    %test if double
    if isa(data_now,'double')
        %test size
        if size(data_now,1) == ypts && size(data_now,2) == xpts
            %assign to prop
            prop.(fields_now{n})=data_now;
        end
    end
end

%convert XSample and YSample to X and Y
try
    prop.x=prop.XSample;
    prop.y=prop.YSample;
    
    prop=fieldclean(prop,'XSample');
    prop=fieldclean(prop,'YSample');
catch
    warning('X Sample and Y Sample Data not converted to x and y Data')
end

%Convert PMap to P
try
    prop.p=prop.PMap;
    prop=fieldclean(prop,'PMap');
    
catch
    warning('No PMap data')
end

%create a grid for x and y indexing
[prop.xg,prop.yg]=meshgrid(1:xpts,1:ypts,1);

%create the orientation data
ori = ...
    rotation('Euler',prop.phi1*degree,prop.PHI*degree,prop.phi2*degree);

%clean the prop
prop=fieldclean(prop,'phi1');
prop=fieldclean(prop,'PHI');
prop=fieldclean(prop,'phi2');

%assign the phase data
if isfield(prop,'Phase')
    phaseassign=prop.Phase;
else
    phaseassign=ones(ypts,xpts);
end

%overwrite the phase data
if nargin > 3
    if ~isempty(phasedata)
        phaseassign=phasedata;
    end
end

%create the EBSD container
ebsd = EBSD(ori, phaseassign,cslist,'options',prop);
end

%remove fields from a structure
function prop=fieldclean(prop,fieldn)
if isfield(prop,fieldn)
    prop=rmfield(prop,fieldn);
end
end


