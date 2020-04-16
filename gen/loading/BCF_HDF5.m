function [ MapData,MicroscopeData,PhaseData,EBSPData ] = BCF_HDF5( InputUser )
%BCF_HDF5 Convert to HDF5 if needed
%build to use the Bruker BCF to HDF5 format converter
%will also load a H5 file & set that up if the file exists
%INPUTS
%InputUser = structure with 
%"EBSD_File" = the BCF/H5 file without file ending
%"BCF_folder" = folder with the BCF located within it
%"HDF5_folder" = folder with H5 (HDF5) or where you would like it made
%
%If the BCF to HDF5 conversion needs to happen, this may take some time
%
%OUTPUTS
%
% MapData = Structure containing EBSD map data - given by pattern number
% MicroscopeData = structure with the phase data contained
% PhaseData = strings with phase information
% EBSPData = information on how to read pattern files

%% Versioning
%v1 - TBB 14/04/2017

if isfield(InputUser,'EBSD_file')
    InputUser.EBSD_File=InputUser.EBSD_file;
end

InputUser.HDF5_file=[InputUser.EBSD_File '.h5'];


InputUser.BCF_file=[InputUser.EBSD_File '.bcf'];

if ~isfield(InputUser,'HDF5_folder')
    InputUser.HDF5_folder=InputUser.Folder;
end

InputUser.HDF5FullFile=fullfile(InputUser.HDF5_folder,InputUser.HDF5_file);


if exist(InputUser.HDF5FullFile) == 0
    InputUser.BCFFullFile=fullfile(InputUser.BCF_folder,InputUser.BCF_file);

    if strcmpi(InputUser.BCFFullFile(end-2:end),'bcf') == 0
        InputUser.BCFFullFile=[InputUser.BCFFullFile '.bcf'];
    end
    
    if exist(InputUser.BCFFullFile) == 0
        error(['The input BCF file does not exist' InputUser.BCFFullFile]);
    end
    
    if exist(InputUser.BCF2HDF5_loc) == 0
        error('The converter tool is absent');
    end
    
    %construct the input arguement
    Text_Run=[InputUser.BCF2HDF5_loc ' "' InputUser.BCFFullFile '" "' InputUser.HDF5FullFile '" "' InputUser.EBSD_File '"'];
    
    %run via DOS
    dos(Text_Run)
end

[ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadHDF5( InputUser );
end

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

EBSP.EDXRaw=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Raw'];
EBSP.EDXCor=[HDF5_info.Groups.Name    '/EDX1/Data/Counts Corrected'];

try
MicroscopeData.EDS_AzimutAngle=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/AzimutAngle']);
MicroscopeData.EDS_Energy=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/Energy']);
MicroscopeData.EDS_PrimaryEnergy=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/PrimaryEnergy']);
MicroscopeData.EDS_TiltAngle=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/TiltAngle']);
MicroscopeData.EDS_NumChan=size(MicroscopeData.EDS_Energy,1);
catch
end

% MicroscopeData.EDS_AzimutAngle=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/AzimutAngle']);
% MicroscopeData.EDS_AzimutAngle=h5read(InputUser.HDF5FullFile,[HDF5_info.Groups.Name    '/EDX1/Header/AzimutAngle']);

%set up the EBSP reader

EBSP.PatternFile=[HDF5_info.Groups.Name    '/EBSD/Data/RawPatterns'];
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

