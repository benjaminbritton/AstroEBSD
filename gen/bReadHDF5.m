function [ MapData,MicroscopeData,PhaseData,EBSP ] = bReadHDF5( InputUser )
%BREADHDF5 Builds EBSD data from a HDF5 file

InputUser.HDF5FullFile=fullfile(InputUser.HDF5_folder,InputUser.HDF5_file);

HDF5_info=h5info(InputUser.HDF5FullFile);

%read the map data
MapData.DD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name         '/EBSD/Data/DD']);
MapData.MAD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/MAD']);
MapData.MADPhase=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name   '/EBSD/Data/MADPhase']);
MapData.NIndexedBands=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/NIndexedBands']);
MapData.PCX=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PCX']);
MapData.PCY=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PCY']);
MapData.PHI=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PHI']);
MapData.phi2=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name       '/EBSD/Data/phi2']);
MapData.phi1=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name       '/EBSD/Data/phi1']);
MapData.RadonQuality=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name '/EBSD/Data/RadonQuality']);
MapData.XBeam=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name      '/EBSD/Data/X BEAM'])+1;
MapData.XSample=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Data/X SAMPLE']);
MapData.YBeam=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name      '/EBSD/Data/Y BEAM'])+1;
MapData.YSample=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Data/Y SAMPLE']);

%read the header info
MicroscopeData.CameraTilt=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/CameraTilt']);
MicroscopeData.KV=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/KV']);
MicroscopeData.MADMax=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MADMax']);
MicroscopeData.Magnification=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Magnification']);
MicroscopeData.MapStepFactor=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MapStepFactor']);
MicroscopeData.MaxRadonBandCount=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MaxRadonBandCount']);
MicroscopeData.NCOLS=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NCOLS']);
MicroscopeData.NROWS=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NROWS']);
MicroscopeData.NPoints=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NPoints']);
MicroscopeData.PatternHeight=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/PatternHeight']);
MicroscopeData.PatternWidth=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/PatternWidth']);
MicroscopeData.SEMImage=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEM Image']);
MicroscopeData.SEPixelSizeX=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEPixelSizeX']);
MicroscopeData.SEPixelSizeY=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEPixelSizeY']);
MicroscopeData.SampleTilt=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SampleTilt']);
MicroscopeData.WD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/WD']);
MicroscopeData.XSTEP=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/XSTEP']);
MicroscopeData.YSTEP=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/YSTEP']);

MicroscopeData.CoordinateSystems=fliplr(rot90(h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Coordinate Systems/ESPRIT Coordinates']),3));

MicroscopeData.CoordSystems='TRZP'; %Top Right Z Plus - standard for Bruker Coords - is used in EBSD_Map

%set up the EBSP reader

EBSP.PatternFile=[HDF5_info.Groups.Name    '/EBSD/Data/RawPatterns'];

EBSP.EDXRaw=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Raw'];
EBSP.EDXCor=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Corrected'];

EBSP.PW=double(MicroscopeData.PatternWidth);
EBSP.PH=double(MicroscopeData.PatternHeight);
EBSP.HDF5_loc=InputUser.HDF5FullFile;

%read the phase info
phase_info=h5info(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases']);

num_phases=size(phase_info.Groups,1);
for n=1:num_phases
    PhaseData(n).Name=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases/' int2str(n) '/Name']);
    PhaseData(n).LatticeConstants=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases/' int2str(n) '/LatticeConstants']);
end
PhaseData(1).num_phases=num_phases;

MicroscopeData.TotalTilt=-((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180; %radians

end
