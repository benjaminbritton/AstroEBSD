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
%v2 - TBB 30/03/2020 - added whole array reading

if exist('pattern_number','var') == 0 %read all the patterns
    patterninfo=h5info(EBSPData.HDF5_loc,EBSPData.PatternFile);
    EBSDPat_t=h5read(EBSPData.HDF5_loc,EBSPData.PatternFile,[1 1 1],[EBSPData.PW EBSPData.PH patterninfo.Dataspace.Size(3)]);
    
    %convert into the array we like to work with
    %swap rows and columns
    EBSDPat=permute(EBSDPat_t,[2,1,3]);
   
    %flip the pattern UD
    EBSDPat=flipud(EBSDPat);
    
    %convert to double
    EBSDPat=double(EBSDPat);
else
    if isnumeric(pattern_number) %check its numeric
        if numel(pattern_number) == 1
            EBSDPat=flipud(double(h5read(EBSPData.HDF5_loc,EBSPData.PatternFile,[1 1 pattern_number],[EBSPData.PW EBSPData.PH 1]))');
        end
        if numel(pattern_number) == 2
           EBSDPat_t=h5read(EBSPData.HDF5_loc,EBSPData.PatternFile,[1 1 pattern_number(1)],[EBSPData.PW EBSPData.PH pattern_number(2)]); 
           EBSDPat=permute(EBSDPat_t,[2,1,3]);
           EBSDPat=flipud(EBSDPat);
           EBSDPat=double(EBSDPat);
        end

    else
        disp('The patterns cannot be read');
    end
end

end

