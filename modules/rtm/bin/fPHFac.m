function [r_fac] = fPHFac(EBSP_geometry,Settings_Cor,screen_int,RTM_setup,SettingsXCF)
%fPHFac - calculate the PH factor correction for a background correction

nr=20;

r_facn=zeros(nr,1);
for n=1:nr
    g_in=rotation.byAxisAngle(vector3d(rand(3,1)),rand(1));
    
    [ EBSP_sim ] = EBSP_gen( EBSP_geometry,g_in.matrix,screen_int); %generate the EBSP for this iteration
    [EBSP_sim_bg] = EBSP_BGCor( EBSP_sim,Settings_Cor );
    
    [Pat_sim_r,XCF_data_fill] = refine_prep(EBSP_sim,SettingsXCF,RTM_setup);
    [Pat_sim_bg_r] = refine_prep(EBSP_sim_bg,SettingsXCF,RTM_setup);
    RegOut_r_rbg= fReg( Pat_sim_r.FFT,Pat_sim_bg_r.FFT,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    r_facn(n)=RegOut_r_rbg(4);
end
r_fac=mean(r_facn(n));

end

