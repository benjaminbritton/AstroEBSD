function [ EBSDPat ] = bReadEBSP(EBSPData,pattern_number)
%BREADEBSP Read an EBSP from a HDF5 file
%EBSPData has:
%
% PatternFile: '/Si_1000x_10um_5_4s/EBSD/Data/RawPatterns'
%              PW: 1600
%              PH: 1152
%        HDF5_loc: 'C:\Users\Benjamin\Documents\Writing\LargeEBSD\Algorithms\Bruker_HDF5\Si_1000x_10um_5_4s.h5'
%
%
% This function has to read the EBSD in the correct coordinate system, or
% everything else will go wrong...

%% Versioning
%v1 - TBB 14/04/2017

EBSDPat=flipud(double(h5read(EBSPData.HDF5_loc,EBSPData.PatternFile,[1 1 pattern_number],[EBSPData.PW EBSPData.PH 1]))');

end

