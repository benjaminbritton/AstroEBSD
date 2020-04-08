function bWriteSpectrum(EBSPData,pattern_number,spectrum_raw,spectrum_cor,channum,Datatype)

%check to see if the array works and needs making
if EBSPData.made == 0
    h5create(EBSPData.HDF5_loc,EBSPData.EDXRaw,[channum,EBSPData.numpats],'Datatype',Datatype);
    h5create(EBSPData.HDF5_loc,EBSPData.EDXCor,[channum,EBSPData.numpats],'Datatype',Datatype);
end

%write the pattern to file
    h5write(EBSPData.HDF5_loc,EBSPData.EDXRaw,spectrum_raw,[1, pattern_number],[channum,1]);
    h5write(EBSPData.HDF5_loc,EBSPData.EDXCor,spectrum_cor,[1, pattern_number],[channum,1]);
end
