function [MapInfo,PatternData] = PC_GridSearch(MapInfo_S,Refine,RTM_setup,Settings_Cor,SettingsXCF,screen_int,RTM_info,Detector_tilt)
%PC_GRIDSEARCH Grid searchers for the pattern centre 
% using RTM gradient descent

x_sel=floor(linspace(2,MapInfo_S.Data_InputMap.xpts-1,Refine.xpc));
y_sel=floor(linspace(2,MapInfo_S.Data_InputMap.ypts-1,Refine.ypc));
[pc_xgrid,pc_ygrid]=meshgrid(x_sel,y_sel);

locations(:,1)=pc_xgrid(:);
locations(:,2)=pc_ygrid(:);

for i=1:numel(pc_xgrid)
    PatternData.PC_start(i,:)=[MapInfo_S.Data_InputMap.PCX(locations(i,2),locations(i,1)),MapInfo_S.Data_InputMap.PCY(locations(i,2),locations(i,1)),MapInfo_S.Data_InputMap.DD(locations(i,2),locations(i,1))]; %initial value for PC
    PatternData.Eulers(i,:)=[MapInfo_S.Data_InputMap.phi1(locations(i,2),locations(i,1))*degree,MapInfo_S.Data_InputMap.PHI(locations(i,2),locations(i,1))*degree,MapInfo_S.Data_InputMap.phi2(locations(i,2),locations(i,1))*degree]; %initial value for Eulers
    %Pattern Number
    PatternData.P(i)=MapInfo_S.Data_InputMap.PMap(locations(i,2),locations(i,1));
    PatternData.Phase(i)=MapInfo_S.Data_InputMap.MADPhase(locations(i,2),locations(i,1));
    if PatternData.Phase(i) == 0
        PatternData.Phase(i) = 1; %overwrite to be the primary phase - could adapt to do SO(3) searching
    end
end

PatternInfo.ScreenWidth=RTM_setup.screensize;
PatternInfo.ScreenHeight=RTM_setup.screensize;

[PatternData.Best_PC,PatternData.Best_PH,PatternData.Best_Ori] = fPCMatchSearch(PatternData,MapInfo_S,Refine,RTM_info,RTM_setup,Settings_Cor,SettingsXCF,screen_int,Detector_tilt,PatternInfo);

% Fit the PC to a plane for the beam positions and create the equivalent map array
[MapInfo] = fPCReCalc(PatternData.Best_PC,PatternData,MapInfo_S);

end

