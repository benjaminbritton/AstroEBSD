clear; home; close all

%set where the data is
InputUser.HDF5_folder='E:\Ben\'; %add trailing '\'
InputUser.HDF5_file= 'Demo_Ben.h5'; %add ending '.h5'

%create the crystal structure list
% cs={'notIndexed',}; %crystal structure list
cs=loadCIF('Fe-Iron-alpha');
% cs={'notIndexed',cs_phase};

[MapInfo_S.MapData,MicroscopeData,PhaseData,MapInfo_S.EBSPData ]=bReadHDF5( InputUser );
[MapInfo_S.Data_InputMap] = EBSD_Map(MapInfo_S.MapData,MicroscopeData);

MTEX_settings.grainangle=5; %grain boundary in degrees
MTEX_settings.xpref='west';
MTEX_settings.ypref='outOfPlane';

%% Create the ebsd container and grains container
[ebsd] = h5_to_MTEX(MapInfo_S.Data_InputMap,MTEX_settings,{'notIndexed',cs});
% calc the grains
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',MTEX_settings.grainangle*degree);

%% Plot both
phase_to_plot=ebsd(1).mineral;

figure;
colorKey = ipfHSVKey(cs(ebsd(1).phase));
colorKey.inversePoleFigureDirection = xvector;
plot(ebsd(phase_to_plot),colorKey.orientation2color(ebsd(phase_to_plot).orientations))
hold on; plot(grains(phase_to_plot).boundary,'color','k','linewidth',1.5); hold off