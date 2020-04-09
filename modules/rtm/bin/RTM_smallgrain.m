function [MapInfo_s] = RTM_smallgrain(grains_small,grains,ebsd,MapInfo_b,Detector_tilt,screen_int_all,RTM_info_all,RTM_setup,Settings_Cor,SettingsXCF)
%RTM_SMALLGRAIN Summary of this function goes here
%   Detailed explanation goes here
MapInfo_s=MapInfo_b;
PatternInfo.ScreenWidth=RTM_setup.screensize;
PatternInfo.ScreenHeight=RTM_setup.screensize;

phi1=MapInfo_s.Data_InputMap.phi1*pi/180;
phi2=MapInfo_s.Data_InputMap.phi2*pi/180;
PHI=MapInfo_s.Data_InputMap.PHI*pi/180;
PH=MapInfo_s.Data_InputMap.PH;

% phi1_s=zeros(size(grains_small,1),1);
% phi2_s=zeros(size(grains_small,1),1);
% PHI_s=zeros(size(grains_small,1),1);
% PH_s=zeros(size(grains_small,1),1);
% xg_s=zeros(size(grains_small,1),1);
% yg_s=zeros(size(grains_small,1),1);

for g=1:size(grains_small,1)
    %extract the neighbour grains
    [neighbour_num,neighbour_gnum]=grains_small(g).neighbors;
    
    %extract properties for these points
    g_p=ebsd(grains_small(g)).prop.p;
    num_g=numel(g_p);
    
    g_pcx=ebsd(grains_small(g)).prop.PCX;
    g_pcy=ebsd(grains_small(g)).prop.PCY;
    g_dd=ebsd(grains_small(g)).prop.DD;
    g_phase=ebsd(grains_small(g)).phase ;
    
    % need to sort this for an orientation rather than the Euler angles
    if g_phase == 1
        g_ori=ebsd(grains_small(g)).orientations.matrix;
    else
        g_ori=repmat(eye(3),1,1,num_g);
    end
    g_ori_det=conv_G_to_Det(g_ori,Detector_tilt);
    
    g_xg=ebsd(grains_small(g)).prop.xg;
    g_yg=ebsd(grains_small(g)).prop.yg;
    
  
    for n=1:num_g %count inside this region
        
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,[g_pcx(n),g_pcy(n),g_dd(n)]);
        [ Pat_Ref ] = bReadEBSP(MapInfo_b.EBSPData,g_p(n));
        %correct this pattern
        [ Pat_Ref_BG] = EBSP_BGCor( Pat_Ref,Settings_Cor );
        [Pat_Ref_r] = refine_prep(Pat_Ref_BG,SettingsXCF,RTM_setup);
        
        trial_phase=grains(neighbour_gnum(:,2)).phase;
        
        trial_ori=grains(neighbour_gnum(trial_phase == 1,2)).meanOrientation;
        trial_phase=trial_phase(trial_phase==1);
        
        num_trial=size(trial_ori,1);
        
        G_F=zeros(3,3,num_trial+1);
        regout_F=zeros(4,num_trial+1);
        
        for t=1:num_trial
            [G_det] = conv_G_to_Det(trial_ori(t).matrix,Detector_tilt);
            
            %needs to be the correct orientation
            [G_F(:,:,t),regout_F(:,t)] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,G_det,SettingsXCF,screen_int_all(trial_phase(t)),RTM_info_all(trial_phase(t)).isHex,RTM_setup);
        end
        
        
             [G_F(:,:,t+1),regout_F(:,t+1)] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,g_ori_det(:,:,n),SettingsXCF,screen_int_all(trial_phase(t)),RTM_info_all(trial_phase(t)).isHex,RTM_setup);
       
        
        %extract max value
        [r_sel,v]=max(regout_F(4,:));
        
        p_ori_sam_refine=conv_G_to_Samp(G_F(:,:,v),Detector_tilt);
        ori_p_sample=rotation('matrix',p_ori_sam_refine);
        
        phi1(g_yg(n),g_xg(n))=ori_p_sample.phi1;
        phi2(g_yg(n),g_xg(n))=ori_p_sample.Phi;
        phi2(g_yg(n),g_xg(n))=ori_p_sample.phi2;
        PH(g_yg(n),g_xg(n))=r_sel;
    end

end

MapInfo_s.Data_InputMap.phi1=phi1*180/pi;
MapInfo_s.Data_InputMap.phi2=phi2*180/pi;
MapInfo_s.Data_InputMap.PHI=PHI*180/pi;
MapInfo_s.Data_InputMap.PH=PH*180/pi;
end

