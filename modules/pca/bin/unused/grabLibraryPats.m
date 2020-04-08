function [pats] = grabLibraryPats(InputUser,RTM,MapInfo,EBSD_Info,Settings_Cor,MicroscopeData)

time1=clock;
screensize=RTM.screensize;
Sampling_Freq=RTM.Sampling_Freq;
XCF_type=RTM.XCF_type;
iterations=RTM.iterations;
LPTsize=RTM.LPTsize;
Bin_loc=RTM.Bin_loc;
Phase_Folder=RTM.Phase_Folder;

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM_PCA( InputUser.Phase_Input,Phase_Folder, Bin_loc );
cs=loadCIF(RTM_info.cif_file);
Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( screensize, LPTsize );
[screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);
[Data_InputMap_Start] = MapInfo.Data_InputMap;
[Data_InputMap]=PC_square(EBSD_Info.EBSPData,Data_InputMap_Start,Settings_Cor);
Detector_tilt = Rx(MicroscopeData.TotalTilt);
PC_X = mean(Data_InputMap.PCX(:));
PC_Y = mean(Data_InputMap.PCY(:));
PC_Z = mean(Data_InputMap.DD(:));

PC_av=[PC_X,PC_Y,PC_Z];
PatternInfo.ScreenWidth=screensize;
PatternInfo.ScreenHeight=screensize;
[ EBSP_av ] = EBSP_Gnom( PatternInfo,PC_av );
[ library_G ] = SO3_rotmat_gen( cs,Sampling_Freq);
[ template_library, pats ] = Library_Gen_PCA(EBSP_av,screen_int,RTM_info.isHex,library_G,XCF_type,SettingsXCF); %Need to update to get rid of XCF_type

end