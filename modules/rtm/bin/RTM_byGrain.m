function [MapInfo_new] = RTM_byGrain(ebsd,grains_big,MapInfo,Settings_Cor,SettingsXCF,RTM_setup,screen_int,RTM_info,t1,Detector_tilt)
% RTM_BYGRAIN use RTM to update the orienations of a grain
% Does not do a SO(3) search
MapInfo_new=MapInfo;

PatternInfo.ScreenWidth=RTM_setup.screensize;
PatternInfo.ScreenHeight=RTM_setup.screensize;

Refine_phi1=MapInfo.Data_InputMap.phi1*degree;
Refine_PHI=MapInfo.Data_InputMap.PHI*degree;
Refine_phi2=MapInfo.Data_InputMap.phi2*degree;
Refine_PH=zeros(size(Refine_phi2));
Refine_rep=zeros(size(Refine_phi2));

pTime(['Starting Big Grain Refinement for ' int2str(size(grains_big,1)) ' grains.'],t1);

for g=1:size(grains_big,1)
    pTime(['Big Grain Refinement, grain = ' int2str(g)],t1);
    %read all the patterns
    
    xpts=ebsd(grains_big(g)).prop.xg;
    ypts=ebsd(grains_big(g)).prop.yg;
    num_p=numel(xpts);
    
    %    ebsd(grains_big).prop.
    [~,ebspBG_tile] = bBlockReadEBSP(xpts,ypts,MapInfo,Settings_Cor);
    
    ebsd_grain=ebsd(grains_big(g));
    
    p_pcx=ebsd_grain.prop.PCX;
    p_pcy=ebsd_grain.prop.PCY;
    p_dd=ebsd_grain.prop.DD;
    
    p_id=ebsd_grain.id;
    
    p_ori=ebsd_grain.orientations.matrix;
    p_ori_det=conv_G_to_Det(p_ori,Detector_tilt);
    
    p_ori_det_refine=zeros(size(p_ori));
    p_PH_out=zeros(num_p,1);
    
    
    
    
    %refine the points within the grain
    
    parfor p=1:num_p
        %prepare the pattern to refine
        [Pat_Ref_r] = refine_prep(ebspBG_tile(:,:,p),SettingsXCF,RTM_setup);
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,[p_pcx(p),p_pcy(p),p_dd(p)]);
        
        %refine
        [p_ori_det_refine(:,:,p),RegOut] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,p_ori_det(:,:,p),SettingsXCF,screen_int,RTM_info.isHex,RTM_setup);
        p_PH_out(p)=RegOut(4);
    end
    p_ori_sam_refine=conv_G_to_Samp(p_ori_det_refine,Detector_tilt);
    ori_p_sample=rotation('matrix',p_ori_sam_refine);
        
    for p=1:num_p
        Refine_phi1(ypts(p),xpts(p))=ori_p_sample(p).phi1;
        Refine_PHI(ypts(p),xpts(p))=ori_p_sample(p).Phi;
        Refine_phi2(ypts(p),xpts(p))=ori_p_sample(p).phi2;
        Refine_PH(ypts(p),xpts(p))=p_PH_out(p);
    end
    
    %search the edge point list
    edge_ID=grains_big(g).boundary.ebsdId;
    has_zero1=find(edge_ID(:,1)==0);
    edge_ID(has_zero1,:)=[];
    
    has_zero2=find(edge_ID(:,2)==0);
    edge_ID(has_zero2,:)=[];
    
    %search the list for points that belong to the grain 
    edge_ID_all=unique(edge_ID);
    in_grain=ismember(edge_ID_all,p_id);
    edge_id_in_grain=edge_ID_all(in_grain);
   
    num_edges=numel(edge_id_in_grain);
    
    xpt=zeros(num_edges,1);
    ypt=zeros(num_edges,1);
    phi1=zeros(num_edges,1);
    PHI=zeros(num_edges,1);
    phi2=zeros(num_edges,1);
    PH=zeros(num_edges,1);
        
    for edge_n=1:num_edges
        edge_lh=find(edge_ID(:,1)==edge_id_in_grain(edge_n));
        edge_rh=find(edge_ID(:,2)==edge_id_in_grain(edge_n));
        edge_neighbour=[edge_ID(edge_lh,2);edge_ID(edge_rh,1)];
        test_set=[edge_neighbour;edge_id_in_grain(edge_n)];
        
        %cheat to only use the single phase
        phase_list=ebsd(test_set).phase;
        test_set(phase_list==0)=[];
        phase_list=ebsd(test_set).phase;
        phase_block=mode(phase_list);
        test_set(phase_list~=phase_block)=[];
        
        
        %extract the orientations and phases of these tests
        ori_test_sample=ebsd(test_set).orientations.matrix;
        phase_test=ebsd(test_set).prop.Phase;
        phase_test(phase_test==0) = 1; %reset this phase if mis indexed
        
        %need the reference pattern for this test
        p=find(p_id==edge_id_in_grain(edge_n));
        
        [Pat_Ref_r] = refine_prep(ebspBG_tile(:,:,p),SettingsXCF,RTM_setup);
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,[p_pcx(p),p_pcy(p),p_dd(p)]);
        
        %refine based upon the edge orientations
        num_test=numel(test_set);
        ori_test_det=conv_G_to_Det(ori_test_sample,Detector_tilt);
        
        p_ori_det_refine=zeros(3,3,num_test);
        p_PH_out=zeros(num_test,1);
        for t=1:num_test
            [p_ori_det_refine(:,:,t),RegOut] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,ori_test_det(:,:,t),SettingsXCF,screen_int(phase_test(t)),RTM_info(phase_test(t)).isHex,RTM_setup);
            p_PH_out(t)=RegOut(4);
        end
        
        [max_ph,v]=max(p_PH_out);
        
        p_ori_sam_refine=conv_G_to_Samp(p_ori_det_refine(:,:,v),Detector_tilt);
        ori_p_sample=rotation('matrix',p_ori_sam_refine);
    
        xpt(edge_n)=ebsd(edge_id_in_grain(edge_n)).xg;
        ypt(edge_n)=ebsd(edge_id_in_grain(edge_n)).yg;
        phi1(edge_n)=ori_p_sample.phi1;
        PHI(edge_n)=ori_p_sample.Phi;
        phi2(edge_n)=ori_p_sample.phi2;
        PH(edge_n)=max_ph;
        
    end
    
    for edge_n=1:num_edges
        Refine_phi1(ypt(edge_n),xpt(edge_n))=phi1(edge_n);
        Refine_PHI(ypt(edge_n),xpt(edge_n))=PHI(edge_n);
        Refine_phi2(ypt(edge_n),xpt(edge_n))=phi2(edge_n);
        Refine_PH(ypt(edge_n),xpt(edge_n))=PH(edge_n);
    end
end

MapInfo_new.Data_InputMap.phi1=Refine_phi1*180/pi;
MapInfo_new.Data_InputMap.PHI=Refine_PHI*180/pi;
MapInfo_new.Data_InputMap.phi2=Refine_phi2*180/pi;
MapInfo_new.Data_InputMap.PH=Refine_PH;
MapInfo_new.Data_InputMap.replaced=Refine_PH;

end

