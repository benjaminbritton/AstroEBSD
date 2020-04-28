function [Best_PC,Best_PH,Best_Ori] = fPCMatchSearch(PatternData,MapInfo,Refine,RTM_info,RTM_setup,Settings_Cor,SettingsXCF,screen_int,Detector_tilt,PatternInfo)
%FPCMATCHSEARCH Take points from a map and find the PC
%               and orientation for these


num_pts=numel(PatternData.P);

Best_PC=zeros(num_pts,3); %the pattern centre for the best fitting data
Best_PH=zeros(num_pts,1); %the quality for the best fitting data
Best_Ori=zeros(3,3,num_pts); %the rotation matrix, in detector
phases=PatternData.Phase;

%note this could be sped up if needed when more phases are being used, as this could be expensive
%screen_int and RTM_info could be sliced better
%likley important to cluster by phase and then parfor on each phase

for n=1:num_pts
    %n
    num_phase=phases(n);
    %extract the pattern from the data
    [ Pat_Ref ] = bReadEBSP(MapInfo.EBSPData,PatternData.P(n));
    %correct this pattern
    [ Pat_Ref_BG] = EBSP_BGCor( Pat_Ref,Settings_Cor );
    [Pat_Ref_r,~] = refine_prep(Pat_Ref_BG,SettingsXCF,RTM_setup);
    
    %refine the orientation
    GMat_test=conv_EA_to_G(PatternData.Eulers(n,:));
    rotmat=GMat_test*Detector_tilt;
    
    if isfield(Refine,'PC_start')
        PCStart=Refine.PC_start;
    else
        PCStart=PatternData.PC_start(n,:);
    end
       
    
    [ EBSD_geom ] = EBSP_Gnom( PatternInfo,PCStart);
    
    %iterate a better orientation
    [RotMat_R] = refine5(Pat_Ref_r,EBSD_geom,PCStart,rotmat,SettingsXCF,screen_int(num_phase),RTM_info(num_phase).isHex,RTM_setup);
    
    %refine the PC for this pattern
    Best_PC(n,:) = fPCRefine(PatternData.PC_start(n,:),Pat_Ref_r,RotMat_R,PatternInfo,Refine,SettingsXCF,screen_int(num_phase),RTM_info(num_phase),RTM_setup);
    
    %update the geometry
    [ EBSD_geom_B ] = EBSP_Gnom( PatternInfo,Best_PC(n,:) );
    %Refine the Orientation
    [Best_Ori(:,:,n),regout] = refine5(Pat_Ref_r,EBSD_geom_B,PatternData.PC_start(n,:),RotMat_R,SettingsXCF,screen_int(num_phase),RTM_info(num_phase).isHex,RTM_setup);
    
    Best_PH(n)=regout(4);
    
    if Refine.debug == 1 %if you need to debug
        Settings_Cor.tkd_onaxis=0; %turns off the resize and centre crop for TKD
        
        %Redrawn from Bruker data
        [ Pat_Template ] = EBSP_gen( EBSD_geom,rotmat,screen_int(num_phase),RTM_info(num_phase).isHex ); %generate the EBSP for this iteration
        [ Pat_Template ] = EBSP_BGCor( Pat_Template,Settings_Cor );
        
        %Refined orientation
        [ Pat_Template_R ] = EBSP_gen( EBSD_geom,RotMat_R,screen_int(num_phase),RTM_info(num_phase).isHex ); %generate the EBSP for this iteration
        [ Pat_Template_R ] = EBSP_BGCor( Pat_Template_R,Settings_Cor ); %correct with the same BG operation
        
        %Better pattern centre
        [ Pat_Template_B ] = EBSP_gen( EBSD_geom_B,Best_Ori(:,:,n),screen_int(num_phase),RTM_info(num_phase).isHex ); %generate the EBSP for this iteration
        [ Pat_Template_B ] = EBSP_BGCor( Pat_Template_B,Settings_Cor );
        
        %Plot the figures
        figure;
        subplot(2,2,1); pPattern(Pat_Ref_BG,EBSD_geom); title('Inp. - Exp.')
        subplot(2,2,2); pPattern(Pat_Template,EBSD_geom); title('Inp. - Sim. BE') %Bruker Euler Angles
        subplot(2,2,3); pPattern(Pat_Template_R,EBSD_geom); title('Out. - Ori. Ref.') %Refined Orientation
        subplot(2,2,4); pPattern(Pat_Template_B,EBSD_geom_B); title('Out. - PC & Ori. Ref.') %Refined Orientation
        
    end
    
end

end

