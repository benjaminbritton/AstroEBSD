function [EBSP,XCF_data_fill] = refine_prep(EBSP_cor,SettingsXCF,RTM_setup)
%REFINE_PREP Prepare a pattern to be refined

screensize=RTM_setup.screensize;
rmin=10;
rmax=screensize/sqrt(2);

EBSP_cor=EBSP_cor-mean(EBSP_cor(:));
    EBSP_cor=EBSP_cor./std(EBSP_cor(:));
% LPTsize=RTM_setup.LPTsize;
[EBSP.FFT,XCF_data_fill]  =fROIEx2(EBSP_cor,SettingsXCF);
EBSP.logp = logsample(EBSP_cor, rmin, rmax, screensize/2, screensize/2, RTM_setup.LPTsize, RTM_setup.LPTsize); %Transform the reference image into LPT space, logsample is in logsample
   
end

