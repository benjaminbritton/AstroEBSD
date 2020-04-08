function bWriteEBSP(EBSPData,pattern_number,EBSP,Datatype)
%BWRITEEBSP Write an EBSP to file

%check to see if the array works and needs making
if EBSPData.made == 0
    h5create(EBSPData.HDF5_loc,EBSPData.PatternFile,[EBSPData.PW EBSPData.PW EBSPData.numpats],'Datatype',Datatype);
end

%write the pattern to file
    h5write(EBSPData.HDF5_loc,EBSPData.PatternFile,transpose(flipud(EBSP)),[1 1 pattern_number],[EBSPData.PW EBSPData.PH 1]);
end

