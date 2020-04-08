function [ Settings_Cor ] = EBSP_StaticBG( Settings_Cor,MicroscopeData,EBSPData)
%EBSP_STATICBG Calculate a static background if needed
%use as [ Settings_Cor ] = EBSP_StaticBG( Settings_Cor,MicroscopeData,EBSPData,t1 )
%
% INPUTS
% Settings_Cor = Background Correction settings, see help for EBSP_BGCor
% MicroscopeData = Structure with information about the scan
%Note this code violates the rule
%that the input should not have the same variable as the output
%(Sorry)

%% Versioning
%v1 - TBB 14/04/2017

%%
if Settings_Cor.RealBG == 1
    if ~isfield(Settings_Cor,'EBSP_bgnum')
        Settings_Cor.EBSP_bgnum=30;
    end
    [EBSP_bg] = EBSP_BGGen( MicroscopeData,EBSPData,Settings_Cor.EBSP_bgnum );
    Settings_Cor.EBSP_bg=EBSP_bg;
end
end

function [ EBSP_bg ] = EBSP_BGGen( MicroscopeData,EBSPData,num_bg )
%EBSP_BGGEN Summary of this function goes here
%   Detailed explanation goes here

max_pats=double(MicroscopeData.NPoints);

patterns_ok=randi(max_pats,num_bg);

EBSP_stack=zeros(MicroscopeData.PatternHeight,MicroscopeData.PatternWidth,num_bg);
BG_Settings.hotpixel=1;
BG_Settings.hot_thresh=1000;
for n=1:num_bg
    EBSP= bReadEBSP(EBSPData,patterns_ok(n));
   [ EBSP_stack(:,:,n),~ ] = EBSP_BGCor( EBSP,BG_Settings );
end

EBSP_bg=sum(EBSP_stack(:,:,:),3)/num_bg;

end