%% Example creation of a HDF5 file
% Takes an existing file and copies it, downsizing the EBSPs
% Not all the variables need copying 
% As now bReadHDF5 gives a warning when groups of data are missing

location_astro='C:\Users\bbrit\Documents\GitHub\AstroEBSD\';
run(fullfile(location_astro,'start_AstroEBSD.m'));

InputUser.HDF5_folder='C:\Users\bbrit\Documents\EBSD';
InputUser.HDF5_file='Demo_Ben.h5';

%read the existing H5 file
[ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadHDF5( InputUser );

%% Write the new one

OutputUser=InputUser;

%new variable names
OutputUser.HDF5_file=[OutputUser.HDF5_file(1:end-3) '_16bin.h5'];
%we assume that this is an h5 file, hdf5 will cause issues with the (end-3)
OutputUser.DataName=InputUser.HDF5_file(1:end-3);
OutputUser.HDF5FullFile=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);

dtype='/EBSD/Data/'; %EBSD data location
htype='/EBSD/Header/'; %Header data location (e.g. microscope settings)

%% Map Data
%pattern centre
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'DD',MapData.DD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCX',MapData.PCX); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCY',MapData.PCY); %write the data

%indexing data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MAD',MapData.MAD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MADPhase',double(MapData.MADPhase)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'NIndexedBands',double(MapData.NIndexedBands)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'RadonQuality',MapData.RadonQuality); %write the data

%orientation data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PHI',MapData.PHI); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi2',MapData.phi2); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi1',MapData.phi1); %write the data

%beam data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X BEAM',double(MapData.XBeam)-1); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y BEAM',double(MapData.YBeam)-1); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X SAMPLE',MapData.XSample); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y SAMPLE',MapData.YSample); %write the data

%% Header Data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MADMax',MicroscopeData.MADMax); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MaxRadonBandCount',MicroscopeData.MaxRadonBandCount); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'CameraTilt',MicroscopeData.CameraTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SampleTilt',MicroscopeData.SampleTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'TotalTilt',MicroscopeData.TotalTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'KV',MicroscopeData.KV); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NCOLS',MicroscopeData.NCOLS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NROWS',MicroscopeData.NROWS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'XSTEP',MicroscopeData.XSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'YSTEP',MicroscopeData.YSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NPoints',MicroscopeData.NPoints); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',MicroscopeData.PatternHeight); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',MicroscopeData.PatternWidth); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'Magnification',MicroscopeData.Magnification); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MapStepFactor',MicroscopeData.MapStepFactor); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeX',MicroscopeData.SEPixelSizeX); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeY',MicroscopeData.SEPixelSizeY); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'WD',MicroscopeData.WD); %write the data

%% Pattern data

binning=16;
dtype='uint16';
EBSPData_1=EBSPData;


[ EBSDPat_1 ] = bReadEBSP(EBSPData_1,1);
EBSDPat_re = imresize(EBSDPat_1,1/binning);
switch dtype
    case 'uint16'
        EBSDPat_re=uint16(EBSDPat_re);
end


EBSPData.HDF5_loc=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);
EBSPData.PW=double(size(EBSDPat_re,2));
EBSPData.PH=double(size(EBSDPat_re,1));
EBSPData.numpats=double(EBSPData.numpats);

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',EBSPData.PH);
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',EBSPData.PW);


disp('Pattern conversion started');
for p=1:EBSPData.numpats
    [ EBSDPat_1 ] = bReadEBSP(EBSPData_1,p);
    EBSDPat_re = imresize(EBSDPat_1,1/binning);
    switch dtype
        case 'uint16'
            EBSDPat_re=uint16(EBSDPat_re);
    end

    if p == 1
        EBSPData.made=0;
        bWriteEBSP(EBSPData,p,EBSDPat_re,dtype)
    else
        EBSPData.made=1;
        bWriteEBSP(EBSPData,p,EBSDPat_re,dtype)
    end
    
    %provide some feedback
    if 1000*round(p/1000) == p
        disp(['Pattern ' int2str(p) ' of ' int2str(EBSPData.numpats) ' patterns converted']);
    end
end

%% Read the data back - warnings will tell you if things have failed
[ MapData_new,MicroscopeData_new,PhaseData,EBSPData ] = bReadHDF5( OutputUser );


%check that we can read this
 [ EBSDPat_back ] = bReadEBSP(EBSPData,1);
