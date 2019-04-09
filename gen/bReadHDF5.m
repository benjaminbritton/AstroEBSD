function [ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadHDF5( InputUser )
%BREADHDF5 Builds EBSD data from a HDF5 file

%create structures
MapData=struct;
MicroscopeData=struct;
PhaseData=struct;
EBSPData=struct;

%create the full data info
InputUser.HDF5FullFile=fullfile(InputUser.HDF5_folder,InputUser.HDF5_file);
EBSPData.HDF5_loc=InputUser.HDF5FullFile;

%check if the H5 is present
if exist(InputUser.HDF5FullFile) == 0
    error('H5 file does not exist');
end

%read info about the h5
HDF5_info=h5info(InputUser.HDF5FullFile);

%set up the EBSP & EDX reader - note that the EDX/EBSD data may not exist!
EBSPData.PatternFile=[HDF5_info.Groups.Name    '/EBSD/Data/RawPatterns'];
EBSPData.EDXRaw=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Raw'];
EBSPData.EDXCor=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Corrected'];


%read the map PC data
try
    MapData.DD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name         '/EBSD/Data/DD']);
    MapData.PCX=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PCX']);
    MapData.PCY=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PCY']);
catch
    warning('PC data not loaded');
end

%read the indexed data
try
    MapData.PHI=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/PHI']);
    MapData.phi2=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name       '/EBSD/Data/phi2']);
    MapData.phi1=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name       '/EBSD/Data/phi1']);
   catch
    warning('Orientation data loading');
end 

%read the associated data
try
    MapData.MAD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/MAD']);
    MapData.MADPhase=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name   '/EBSD/Data/MADPhase']);
    MapData.NIndexedBands=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name        '/EBSD/Data/NIndexedBands']);
    MapData.RadonQuality=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name '/EBSD/Data/RadonQuality']);
    MicroscopeData.MaxRadonBandCount=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MaxRadonBandCount']);
    MicroscopeData.MADMax=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MADMax']);
catch
    warning('Indexing of patterns not loaded')
end

try
    MapData.PeakHeight=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name         '/EBSD/Data/PeakHeight']);
catch
    warning('PeakHeight data not loaded');
end

%read the beam data
try
    MapData.XBeam=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name      '/EBSD/Data/X BEAM'])+1; %note that 1 is added
    MapData.XSample=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Data/X SAMPLE']);
    MapData.YBeam=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name      '/EBSD/Data/Y BEAM'])+1; %note that 1 is added
    MapData.YSample=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Data/Y SAMPLE']);
    MapData.PMap=transpose(1:numel(MapData.YSample));
catch
    warning('Beam and sample data not loaded')
end
%read the header info

try
    MicroscopeData.CameraTilt=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/CameraTilt']);
catch
    warning('Camera tilt not loaded');
end

try
    MicroscopeData.SampleTilt=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SampleTilt']);
catch
    warning('Sample tilt not loaded')
end

if isfield(MicroscopeData,'SampleTilt') && isfield(MicroscopeData,'CameraTilt')
    MicroscopeData.TotalTilt=-((90-MicroscopeData.SampleTilt)+MicroscopeData.CameraTilt)*pi/180; %radians
else
    MicroscopeData.TotalTile=NaN;
end

try
    MicroscopeData.KV=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/KV']);
catch
    warning('KV not loaded');
end

try
    MicroscopeData.NCOLS=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NCOLS']);
    MicroscopeData.NROWS=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NROWS']);
    MicroscopeData.XSTEP=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/XSTEP']);
    MicroscopeData.YSTEP=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/YSTEP']);
    MicroscopeData.NPoints=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/NPoints']);
    EBSPData.numpats=MicroscopeData.NPoints;
catch
    warning('Map scanning data not known')
end

try
    MicroscopeData.PatternHeight=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/PatternHeight']);
    MicroscopeData.PatternWidth=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/PatternWidth']);
    EBSPData.PW=double(MicroscopeData.PatternWidth);
    EBSPData.PH=double(MicroscopeData.PatternHeight);
catch
    warning('Pattern information not loaded')
end

try
    MicroscopeData.Magnification=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Magnification']);
    MicroscopeData.MapStepFactor=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/MapStepFactor']);

    MicroscopeData.SEPixelSizeX=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEPixelSizeX']);
    MicroscopeData.SEPixelSizeY=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEPixelSizeY']);
    MicroscopeData.WD=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/WD']);
catch
    warning('Microscope settings data not loaded');
end

try
    MicroscopeData.SEMImage=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/SEM Image']);
catch
    warning('SE Image not loaded');
end

try
    MicroscopeData.CoordinateSystems=fliplr(rot90(h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Coordinate Systems/ESPRIT Coordinates']),3));
catch
    warning('Coordinate systems not loaded')
end
    MicroscopeData.CoordSystems='TRZP'; %Top Right Z Plus - standard for Bruker Coords - is used in EBSD_Map

    

%read the phase info
try
    phase_info=h5info(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases']);
    
    num_phases=size(phase_info.Groups,1);
    for n=1:num_phases
        PhaseData(n).Name=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases/' int2str(n) '/Name']);
        PhaseData(n).LatticeConstants=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EBSD/Header/Phases/' int2str(n) '/LatticeConstants']);
    end
    PhaseData(1).num_phases=num_phases;
catch
    warning('Phase data not loaded')
end


end
