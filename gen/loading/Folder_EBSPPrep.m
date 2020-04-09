function EBSP_One=Folder_EBSPPrep(FolderOut,InputUser,Settings_Cor,Settings_Rad,n,Num_Phase)
%read the pattern
EBSP_One.PatternIn= ReadEBSDFile(FolderOut.patternlist{n},InputUser.PatternFlip);

EBSP_One.PC=FolderOut.PC_out(:,1,n);

[ EBSP_One.PatternCor,EBSP_One.PatternInfo ] = EBSP_BGCor( EBSP_One.PatternIn,Settings_Cor );
% radon convert & Peak ID
[ EBSP_One.Peak_Centre,EBSP_One.Single.Peak_Set_All,EBSP_One.Peak_Set_All,...
    EBSP_One.R_EBSP,EBSP_One.R_Edge,EBSP_One.R_rho,EBSP_One.R_theta ] ...
    = EBSP_RadHunt( EBSP_One.PatternCor,Settings_Rad);
[ EBSP_One.nhat_gnom] = EBSP_NormConv( EBSP_One.Peak_Centre,[EBSP_One.PatternInfo.size],EBSP_One.PC);
            
[ EBSP_One.PatternGeometry ] = EBSP_Gnom( EBSP_One.PatternInfo,EBSP_One.PC );
   
EBSP_One.rotdata{Num_Phase}.detector=FolderOut.rotdata{1,Num_Phase}.detector;