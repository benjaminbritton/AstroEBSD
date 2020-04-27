function [MapInfo2]=fRefineMap(Refine,RTM,MapInfo,MicroscopeData,InputUser,Settings_Cor,t1)

[ SettingsXCF] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );
PatternInfo.ScreenWidth=RTM.screensize;
PatternInfo.ScreenHeight=RTM.screensize;
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

%% use RTM (non SO(3) search) for the PC Searching
PhaseInput=InputUser.Phases{Refine.phase};
pTime(['Set up RTM PC refinement to pattern match against ' PhaseInput],t1);
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( {PhaseInput},RTM.Phase_Folder);
%generate the simulated pattern cube
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Select points and refine

%pTime(['Point selection for PC refinement for ' int2str(Refine.num_pts) ' points'],t1);

[PatternData] = fSelectPoints(MapInfo,PhaseInput,Refine.num_pts);
close all

%tell PatternData what phase to search for for each point
PatternData.Phase=repmat(Refine.phase,[Refine.num_pts,1]);

pTime(['Optimising pattern centre for selected points'],t1);
% Find the PC for these points
[Best_PC,Best_PH,Best_Ori] = fPCMatchSearch(PatternData,MapInfo,Refine,RTM_info,RTM,Settings_Cor,SettingsXCF,screen_int,Detector_tilt,PatternInfo);

% Fit the PC to a plane for the beam positions and create the equivalent map array
[MapInfo2] = fPCReCalc(Best_PC,PatternData,MapInfo);
pTime(['Refinement complete'],t1);
MapInfo2.MapInfo_Original=MapInfo;

%add a flag to say we've refined the map
MapInfo2.Refined=1;
end
