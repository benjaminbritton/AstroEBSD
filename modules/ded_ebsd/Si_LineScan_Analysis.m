%% Si_LineScan_Analysis.m
% Originally by Tianbi Zhang, November 2023
% Edited by T Ben Britton
% This code analyzes the EBSPs from the line scan experiments on the
% Si(001) single crystal sample. It reads background corrected patterns
% (image files) and analyze one by one.

% Requirements:
% Si patterns (see doi:10.5281/zenodo.8342489)
% MTEX toolbox
% MATLAB image processing, curve fitting and statisics and machine learning
% toolboxes

% For fitting the PC values to beam shift, please check the 'static_EBSD'
% folder.

%% Start of the script
close all;
clear;
home;

%% start astro & mtex - please edit the paths!

InputUser.Astro_loc = 'C:\Users\billy\Documents\GitHub\AstroEBSD';
InputUser.location_mtex='E:\MATLAB\mtex-5.8.0\';
run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m'));
run(fullfile(InputUser.location_mtex,'startup.m'));

pattern_loc = 'C:\Users\billy\OneDrive - UBC\PhD\TKD\20230804_Si_EBSD\Vertscan';

%% Enter pattern information and set up indexing
eangs=[176.81 19.31 225.90]*pi/180; %demo
Settings_PCin.start=[0.4595856 0.5319605 0.7521288];% demo - vertnew
% Settings_PCin.start=[0.4723365 0.5362904 0.7578808];% demo - horinew
Settings_PCin.range = [0.0002 0.0002 0.0002];

MicroscopeData.TotalTilt=0;
%Set up the radon transform peak finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak huntfor the Hough map
Settings_Rad.max_peaks=12; %max number of peaks to return
Settings_Rad.num_peak=30; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=6; %search size in theta steps
Settings_Rad.rho_search_per=0.12; %radon search in fractions
Settings_Rad.min_peak_width=0.002; %min rseperation of the peak width, in pixels

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

% gmatrix = RTM.Rz(eangs(3)) * RTM.Rx(eangs(2)) * RTM.Rz(eangs(1));

R_x=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
rot_det=R_x(MicroscopeData.TotalTilt);

%% Pattern I/O
num_pat_ini = 1;
num_pat_end = 1;
pattern_root = 'Spot';
num_pats = num_pat_end - num_pat_ini + 1;

eangs_list = zeros(3,num_pats);
PC_list = zeros(num_pats,3);
PCx_est = Settings_PCin.start(1) - 10 / 252 / 55 * (num_pat_ini:num_pat_end)';
PCz_est = Settings_PCin.start(3) + 0 /252 /55 .* ((num_pat_ini-1:num_pat_end-1)');
PCy_est = Settings_PCin.start(2) - (tan(70 * pi / 180) * 0 / 252 / 55) .* ((num_pat_ini-1:num_pat_end-1)'); 


InputUser.Phase_Folder = fullfile(InputUser.Astro_loc,'phases\phasefiles');
InputUser.Phase_Input  = {'Si_30kv_HR'};
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );
cs = crystalSymmetry('m-3m'); ss = specimenSymmetry('orthorhombic');

Settings_Cor.radius=1;
Settings_Cor.radius_frac=0.5;
Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=20; %low pass filter sigma
Settings_CorX.hotpixel=0; %hot pixel correction
Settings_CorX.hot_thresh=1000; %hot pixel threshold

%set up the PC refinement algorithm
Refine.ss=0.08; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=50;
Refine.reindex=1;
Refine.print=1;
Detector_tilt = eye(3); %RTM.Rx(20 * pi / 180);
% [EBSD_geom_init ] = EBSP_Gnom( RTM,Settings_PCin.start); %you can change PC_in if you want

%% Load the EBSPs
disp('Loading EBSPs into RAM')
pattern_all=zeros(252,252,400); 
for pat=num_pat_ini:num_pat_end
    pattern_name=[pattern_root int2str(pat)];
    pattern_read = hdfpatternread(pattern_loc, pattern_name, true);
    [pattern_all(:,:,pat),~] = EBSP_BGCor( pattern_read,Settings_Cor);
end

%% index in AstroEBSD (Hough)

for i=num_pat_ini:num_pat_end

pattern_name = strcat(pattern_root,int2str(i),'_bgcor.tif');
pattern_file = fullfile(pattern_loc,pattern_name);
pat = double(flipud(imread(pattern_file)));

% Normalise intensities
[EBSP_Fe_PCOut.PatternIn,Settings_Cor ] = EBSP_BGCor( pat,Settings_Cor);
EBSP_Fe_PCOut.PC = Settings_PCin.start;

%Radon transform
[ EBSP_Fe_PCOut.Peak_Centre,EBSP_Fe_PCOut.Single.Peak_Set_All,EBSP_Fe_PCOut.Peak_Set_All,...
    EBSP_Fe_PCOut.R_EBSP,EBSP_Fe_PCOut.R_Edge,EBSP_Fe_PCOut.R_rho,EBSP_Fe_PCOut.R_theta ] ...
    = EBSP_RadHunt( EBSP_Fe_PCOut.PatternIn,Settings_Rad);

% Convert the bands to normal space
[ EBSP_Fe_PCOut.nhat_gnom] = EBSP_NormConv( EBSP_Fe_PCOut.Peak_Centre,size(EBSP_Fe_PCOut.PatternIn),EBSP_Fe_PCOut.PC);

 % Astro_EBSPset(EBSP_One,Settings_Cor,Settings_Rad,pc_pat1,InputUser);

[EBSP_Fe_PCOut] = EBSP_PCSearch(EBSP_Fe_PCOut,Settings_Cor,Settings_Rad,Settings_PCin,Phase_Num,Crystal_LUT,Crystal_UCell);

[ EBSP_Fe_PCOut.nhat_gnom] = EBSP_NormConv( EBSP_Fe_PCOut.Peak_Centre,size(EBSP_Fe_PCOut.PatternIn),EBSP_Fe_PCOut.PC);
[EBSP_Fe_PCOut.rotdata{1},EBSP_Fe_PCOut.banddata{1}]=EBSP_Index(EBSP_Fe_PCOut.nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},rot_det); %#ok<PFBNS>

%generate the geometry
[ EBSP_Fe_PCOut.PatternGeometry ] = EBSP_Gnom( Settings_Cor,EBSP_Fe_PCOut.PC);

EBSP_OneFigure=Plot_SinglePattern(EBSP_Fe_PCOut,Crystal_UCell,Crystal_LUT,1);


eangs_list(:,i-num_pat_ini+1) = EBSP_Fe_PCOut.rotdata{1}.eang;
error_list(i-num_pat_ini+1) = EBSP_Fe_PCOut.rotdata{1}.error;
PC_list(i-num_pat_ini+1,:) = EBSP_Fe_PCOut.PC_out';

% change these lines for vertical scan as they update the PC guess.
Settings_PCin.start(1) = Settings_PCin.start(1) - 10 / 252 / 55;
% Settings_PCin.start(2) = EBSP_Fe_PCOut.PC(2) - tan(70 * pi / 180) * 10 / 252 / 55;
% Settings_PCin.start(3) = EBSP_Fe_PCOut.PC(3) + 10 / 252 / 55;

end

eangs_mtex = eangs_list';
ori = orientation.byEuler(eangs_mtex(:,1),eangs_mtex(:,1),eangs_mtex(:,1),'ZYZ',cs,ss);

%% Plots of MAE and PC shift

figure;
plot(ori,'filled');
figure;
hold on;
subplot(2,1,1);
histogram(error_list, 'EdgeColor','none');
xlabel('MAE (degrees)');
ylabel('Counts');
subplot(2,1,2);
plot(PC_list(:,1)); hold on; plot(PC_list(:,2)); plot(PC_list(:,3));
plot(PCx_est,'.'); plot(PCy_est,'.'); plot(PCz_est,'.'); % for vertscan
% plot(PCx_est,'.'); % for horiscan
legend('PCx', 'PCy', 'PCz', 'PCx exp.', 'PCy exp.','PCz exp.', 'Location','eastoutside');
grid on;
xlabel('Pattern #');
ylabel('PCx, PCy, PCz');

pat_id = num_pat_ini:num_pat_end;
pcx_list=PC_list(:,1);
pcy_list=PC_list(:,2);
pcz_list=PC_list(:,3);
