function EBSP_One=Pattern_fromMap(EBSPposition,PCOut,EBSD_DataInfo,Data_OutputMap,Settings_Cor,Settings_Rad,Crystal_LUT,Settings_LUT,Crystal_UCell,Phase_Num)
%% Index a Single Pattern

%choose a point from a map
EBSP_One.x_co=EBSPposition(1); EBSP_One.y_co=EBSPposition(2);
EBSP_One.pat_num=Data_OutputMap.P_co(EBSP_One.y_co,EBSP_One.x_co); %row,column, i.e. y,x
EBSP_One.PC=[PCOut.Fit_2nd.PCx_map(EBSP_One.y_co,EBSP_One.x_co) PCOut.Fit_2nd.PCy_map(EBSP_One.y_co,EBSP_One.x_co) PCOut.Fit_2nd.PCz_map(EBSP_One.y_co,EBSP_One.x_co)];
% Single.PC_data2(Single.pat_num,:)=Single.PC_data2(Single.pat_num-2,:);

%load the pattern
EBSP_One.PatternIn = bReadEBSP(EBSD_DataInfo,EBSP_One.pat_num);
%bg correct
[ EBSP_One.PatternCor,EBSP_One.PatternInfo ] = EBSP_BGCor( EBSP_One.PatternIn,Settings_Cor );

%generate the geometry
[ EBSP_One.PatternGeometry ] = EBSP_Gnom( EBSP_One.PatternInfo,EBSP_One.PC );

% radon convert & Peak ID
[ EBSP_One.Peak_Centre,EBSP_One.Single.Peak_Set_All,EBSP_One.Peak_Set_All,EBSP_One.R_EBSP,EBSP_One.R_Edge,EBSP_One.R_rho,EBSP_One.R_theta ] = EBSP_RadHunt( EBSP_One.PatternCor,Settings_Rad);

% Convert the bands to normal space
[ EBSP_One.nhat_gnom] = EBSP_NormConv( EBSP_One.Peak_Centre,[EBSP_One.PatternInfo.size],EBSP_One.PC);

rot_det=eye(3);
%index for this phases
[EBSP_One.rotdata,EBSP_One.banddata]=EBSP_Index(EBSP_One.nhat_gnom,Crystal_LUT{Phase_Num},Settings_LUT{Phase_Num}.thresh_trig,Crystal_UCell{Phase_Num},rot_det); %#ok<PFBNS>