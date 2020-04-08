function [EDSData_cor,EDSData_raw ] = bReadEDX(EBSPData,pattern_number,channum)
%BREADEBSP Read an EDX data from
%EBSPData has:
%
%
% This function is hard coded to read 2048 spectrum channels

%% Versioning
%v1 - TBB 14/04/2017

if exist('pattern_number','var') == 0
    eds_info=h5info(EBSPData.HDF5_loc,EBSPData.PatternFile);
            EDSData_raw=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXRaw,[1 1],[channum eds_info.Dataspace.Size(3)]));
            EDSData_cor=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXCor,[1 1],[channum eds_info.Dataspace.Size(3)]));
            
else
    if isnumeric(pattern_number) %check its numeric
        if numel(pattern_number) == 1
            EDSData_raw=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXRaw,[1 pattern_number],[channum 1]));
            EDSData_cor=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXCor,[1 pattern_number],[channum 1]));
        else
            EDSData_raw=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXRaw,[1 pattern_number(1)],[channum pattern_number(2)]));
            EDSData_cor=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXCor,[1 pattern_number(1)],[channum pattern_number(2)]));
        end
    end

end

