function [EDSData_cor,EDSData_raw ] = bReadEDX(EBSPData,pattern_number,channum)
%BREADEBSP Read an EDX data from
%EBSPData has:
%
%
% This function is hard coded to read 2048 spectrum channels

%% Versioning
%v1 - TBB 14/04/2017

EDSData_raw=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXRaw,[1 pattern_number],[channum 1]));
EDSData_cor=double(h5read(EBSPData.HDF5_loc,EBSPData.EDXCor,[1 pattern_number],[channum 1]));


end

