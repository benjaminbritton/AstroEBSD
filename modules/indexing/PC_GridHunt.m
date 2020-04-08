function [PCfun_x,PCfun_y,PCfun_z]=PC_GridHunt(PC_grid,Data_InputMap,EBSD_DataInfo,Settings_CorA,Crystal_LUT,Crystal_UCell,Settings_Rad)
[PC_xgrid,PC_ygrid]=meshgrid(1:PC_grid(2),1:PC_grid(1));
PC_xgrid=floor((Data_InputMap.xpts/(PC_grid(2)+1))*PC_xgrid);
PC_ygrid=floor((Data_InputMap.ypts/(PC_grid(1)+1))*PC_ygrid);

PC_xgrid=PC_xgrid(:);
PC_ygrid=PC_ygrid(:);

num_grid=size(PC_xgrid,1);

PC_GA_options = optimoptions('ga');
PC_GA_options.FunctionTolerance=1E-3;
PC_GA_options.UseParallel=0;
PC_GA_options.MaxGenerations=15;
PC_GA_options.PopulationSize=30;
PC_GA_options.MaxStallGenerations=20;
PC_GA_options.Display='off';
PC_foundp=zeros(3,num_grid,size(Crystal_LUT,2));
PC_errorp=zeros(1,num_grid,size(Crystal_LUT,2));

for num_P=1:size(Crystal_LUT,2)
    PC_found=zeros(3,num_grid);
    PC_error=zeros(1,num_grid);
    parfor g=1:num_grid
        pnum=Data_InputMap.PMap(PC_ygrid(g),PC_xgrid(g));
        EBSP_raw = bReadEBSP(EBSD_DataInfo,pnum);
        
        Settings_PCin_start=[Data_InputMap.PCX(pnum) Data_InputMap.PCY(pnum) Data_InputMap.DD(pnum)];
        Settings_PCin_range=[0.1 0.1 0.1];
        
        [ PatternCor,PatternInfo ] = EBSP_BGCor( EBSP_raw,Settings_CorA );
        % radon convert & Peak ID
        [ EBSP_Onep_Peak_Centre,EBSP_Onep_Single_Peak_Set_All,EBSP_Onep_Peak_Set_All,...
            EBSP_Onep_R_EBSP,EBSP_Onep_R_Edge,EBSP_Onep_R_rho,EBSP_Onep_R_theta ] ...
            = EBSP_RadHunt( PatternCor,Settings_Rad);
        % find the pattern centre
        
        PC_GA_ub=Settings_PCin_start+Settings_PCin_range;
        PC_GA_lb=Settings_PCin_start-Settings_PCin_range;
%         Phase_Num=1;
%         EBSP_Onep_PC_out=zeros(3,Phase_Num);
%         EBSP_Onep_PC_err=zeros(Phase_Num,1);
%         
        FitFunc = @(PC_test) PC_GAOpt( PC_test,EBSP_Onep_Peak_Centre,PatternInfo.size,Crystal_LUT,Crystal_UCell,num_P);
        [EBSP_Onep_PC_out, EBSP_Onep_PC_err] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
        
        PC_found(:,g)=EBSP_Onep_PC_out;
        PC_error(g)=EBSP_Onep_PC_err;
      
    end
    PC_foundp(:,:,num_P)=PC_found;
    PC_errorp(:,:,num_P)=PC_error;
end

%now find the phase for each

[PC_error_f,ix]=min(PC_errorp,[],3);
PC_found_f=0*PC_found;
PC_thresh=0.1;
for g=1:num_grid
    PC_found_f(:,g)=PC_foundp(:,g,ix(g));
    if PC_error_f > PC_thresh
        PC_found_f(:,g)=NaN*PC_found_f(:,g);
    end
end
xyg=[PC_xgrid PC_ygrid];
PCfun_x=robustfit(xyg,PC_found_f(1,:));
PCfun_y=robustfit(xyg,PC_found_f(2,:));
PCfun_z=robustfit(xyg,PC_found_f(3,:));

% %to eval
% PC_X=PCfun_x(1)+Data_InputMap.XBeam_Map*PCfun_x(2)+Data_InputMap.YBeam_Map*PCfun_x(3);
% PC_Y=PCfun_y(1)+Data_InputMap.XBeam_Map*PCfun_y(2)+Data_InputMap.YBeam_Map*PCfun_y(3);
% PC_Z=PCfun_z(1)+Data_InputMap.XBeam_Map*PCfun_z(2)+Data_InputMap.YBeam_Map*PCfun_z(3);
