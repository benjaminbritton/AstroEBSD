%% stereo_reproject.m
% Created by T. Ben Britton in June 2023
% Edited and annotated by Tianbi Zhang

% This script takes a Kikuchi pattern, index it, reproject it onto the
% diffraction sphere and create a stereogram.

% This script accompanies the following article:
% "Multi-exposure diffraction pattern fusion applied to enable wider-angle 
% transmission Kikuchi diffraction with direct electron detectors"
% Tianbi Zhang, T. Ben Britton 

% Requirements: 
% (1) MATLAB toolboxes: image processsing, statistics and machine learning,
% parallel computing
% (2) AstroEBSD package - this script is available as a part of the latest
% AstroEBSSD distribution.
% (3) MTEX toolbox (https://mtex-toolbox.github.io/)

%% Start of the script
close all;
clear;
home;

%% start astro & mtex - please edit the paths
location_astro='C:\Users\benja\OneDrive\Documents\GitHub\AstroEBSD\';
location_mtex='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.4.0';
path1= 'C:\Users\benja\OneDrive\Documents\MATLAB\Al_HDR\Demo';


run(fullfile(location_astro,'start_AstroEBSD.m'));
run(fullfile(location_mtex,'startup.m'));


%% load the pattern + Figure 1


pat1=double(flipud(imread(fullfile(path1,'stitched_flat_gaus.tif')))); % this is the demo TKP
figure;
imagesc(pat1); axis image; axis xy; colormap('gray');

%% Enter pattern information and set up indexing

% Euler angles (Bunge) and PC - they must be very accurate for the
% reprojection.
eangs_pat1=[140.4739,23.8782,246.8339];%demo
pc_pat1=[0.5454366,0.4350097,0.7578597];% demo

% eangs_pat1=[277.5706,40.1078,68.0877];%ap1
% pc_pat1=[0.5308719,0.4267852,0.6428683];% ap1

% eangs_pat1=[351.0974,35.4311,350.7168];%ap2
% pc_pat1=[0.5402837,0.4316793,0.6424058];% ap2

MicroscopeData.TotalTilt=0;
%Set up the radon transform peak finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak huntfor the Hough map
Settings_Rad.max_peaks=12; %max number of peaks to return
Settings_Rad.num_peak=20; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=6; %search size in theta steps
Settings_Rad.rho_search_per=0.2; %radon search in fractions
Settings_Rad.min_peak_width=0.002; %min rseperation of the peak width, in pixels


% Normalise intensities
[EBSP_One.PatternIn,Settings_Cor ] = EBSP_BGCor( pat1,[]);

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%% index in AstroEBSD (Hough)

EBSP_Fe_One.PatternIn=pat1;
EBSP_Fe_One.PC=pc_pat1;

InputUser.Phase_Folder = fullfile(location_astro,'phases\phasefiles');
InputUser.Phase_Input  = {'Al'};
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );

%Radon transform
[ EBSP_Fe_One.Peak_Centre,EBSP_Fe_One.Single.Peak_Set_All,EBSP_Fe_One.Peak_Set_All,...
    EBSP_Fe_One.R_EBSP,EBSP_Fe_One.R_Edge,EBSP_Fe_One.R_rho,EBSP_Fe_One.R_theta ] ...
    = EBSP_RadHunt( EBSP_Fe_One.PatternIn,Settings_Rad);

% Convert the bands to normal space
[ EBSP_Fe_One.nhat_gnom] = EBSP_NormConv( EBSP_Fe_One.Peak_Centre,size(EBSP_Fe_One.PatternIn),EBSP_Fe_One.PC);

R_x=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
rot_det=R_x(MicroscopeData.TotalTilt);

% Index this pattern
[EBSP_Fe_One.rotdata{1},EBSP_Fe_One.banddata{1}]=EBSP_Index(EBSP_Fe_One.nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},rot_det);

%generate the geometry
[ EBSP_Fe_One.PatternGeometry ] = EBSP_Gnom( Settings_Cor,EBSP_Fe_One.PC );

EBSP_OneFigure=Plot_SinglePattern(EBSP_Fe_One,Crystal_UCell,Crystal_LUT,1);

%% Simulate the pattern

eangs=eangs_pat1*pi/180;
PC_simulation=pc_pat1;
pattern_info=struct;
pattern_info.size=size(pat1);
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

% InputUser.Phase_Folder
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},[location_astro 'phases']);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

[EBSD_simulation ] = EBSP_Gnom( pattern_info,PC_simulation); %you can change PC_in if you want
gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
[ Pat_sim_eang ] = EBSP_gen( EBSD_simulation,gmatrix*Detector_tilt,screen_int); %generate the EBSP for this iteration

[Pat_sim_astroind]= EBSP_gen( EBSD_simulation,EBSP_Fe_One.rotdata{1}.detector,screen_int); %generate the EBSP for this iteration

%% refine + Figure 2

RTM.screensize = 254; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM.iterations = 10;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 254; %LPT size used in pixels
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );

[ Pat_Ref_r ] = refine_prep(EBSP_Fe_One.PatternIn,SettingsXCF,RTM);
%perform a refinement
[G_refined,regout_Grefined] = refine5(Pat_Ref_r,EBSD_simulation,EBSD_simulation.PC,gmatrix,SettingsXCF,screen_int,RTM_info.isHex,RTM);
[ Pat_sim_ref ] = EBSP_gen( EBSD_simulation,G_refined,screen_int); %generate the EBSP for this iteration

figure;
subplot(1,2,1);
pPattern(Pat_sim_ref,EBSD_simulation);
subplot(1,2,2);
pPattern(EBSP_Fe_One.PatternIn,EBSD_simulation);

%% Figure 3 - fitted/indexed patterns
figure;

subplot(1,3,1);
I1=pPattern(pat1,EBSD_simulation);
title('Experimental pattern');

subplot(1,3,2);
I1=pPattern(Pat_sim_eang,EBSD_simulation);
title('Simulated pattern - Dynamics');

subplot(1,3,3);
I1=pPattern(Pat_sim_astroind,EBSD_simulation);
title('Simulated pattern - AstroH');

%% Now try to plot a sterogram of the dynamically simulated pattern + Figure 4
num_stereo=2001;
stereo_range=1;
stereo_xl=linspace(-stereo_range,stereo_range,num_stereo);
stereo_yl=linspace(-stereo_range,stereo_range,num_stereo);

[stereo_x,stereo_y]=meshgrid(stereo_xl,stereo_yl);

stereo_r=stereo_x.^2+stereo_y.^2;

stereo_nu=2./(stereo_r+1);

stereo_vx=stereo_nu.*stereo_x;
stereo_vy=stereo_nu.*stereo_y;
stereo_vz=-1+stereo_nu;

stereo_vl=sqrt(stereo_vx.^2+stereo_vy.^2+stereo_vz.^2);
stereo_vx=stereo_vx./stereo_vl;
stereo_vy=stereo_vy./stereo_vl;
stereo_vz=stereo_vz./stereo_vl;


%sample from the sphere
[stereo_pat] = Cube_Sample(stereo_vx(:),stereo_vy(:),stereo_vz(:),screen_int,0);
% figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat,num_stereo,num_stereo));
% axis image; axis xy; colormap('gray');

%% Now take the experimental pattern and generate a look up, so we can do the spherical sampling...
% Figure 5
% rotmat=EBSP_Fe_One.rotdata{1}.detector;
rotmat=gmatrix*Detector_tilt;
[Pat_sim_astroind,r_exp]= EBSP_gen( EBSD_simulation,rotmat,screen_int); %generate the EBSP for this iteration

r_lambda=1./(r_exp(:,3)+1);

%convert from r_exp to the stereogram
P_rx=r_exp(:,1).*r_lambda;
P_ry=r_exp(:,2).*r_lambda;

pat1_plot=pat1;
pat1_plot=pat1_plot/std(pat1(:));
pat1_plot=pat1_plot-mean(pat1_plot(:));

% create a intensity normalized stereogram
stereo_pat_plot=stereo_pat;
stereo_pat_plot=stereo_pat/std(stereo_pat_plot(:));
stereo_pat_plot=stereo_pat_plot-mean(stereo_pat_plot(:));
%
% figure;
% hist(pat1_plot(:),100,'r');
% hold on;
% hist(stereo_pat_plot(:),100,'b');

figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray'); axis off;
hold on;
% here the pattern is plotted as a scatter plot whose marker color scale
% corresponds to the grayscale value. This is still OK at low mag.
% scatter(P_rx,P_ry,1,pat1_plot(:)*0.8);

%% use the interpolant to generate the pattern
% create the interpolant
stero_interp_exp=scatteredInterpolant(P_rx(:),P_ry(:),pat1_plot(:),'natural','none');

%read the intepolant
stereo_remap=stero_interp_exp(stereo_x(:),stereo_y(:));

stereo_remap_2d=reshape(stereo_remap,num_stereo,num_stereo);
stereo_mask_2d=isnan(stereo_remap_2d);

figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray');
hold on;
i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2d);
axis image; axis xy; colormap('gray');
i1.AlphaData=1-stereo_mask_2d;

% plot the outline of the detector (experimental pattern)
line_top=[EBSD_simulation.xpts_screen(:,1),EBSD_simulation.ypts_screen(:,1)];
line_bottom=flipud([EBSD_simulation.xpts_screen(:,end),EBSD_simulation.ypts_screen(:,end)]);
line_right=[EBSD_simulation.xpts_screen(end,:);EBSD_simulation.ypts_screen(end,:)]';
line_left=flipud([EBSD_simulation.xpts_screen(1,:);EBSD_simulation.ypts_screen(1,:)]');

line_allx=[line_top(:,1);line_right(:,1);line_bottom(:,1);line_left(:,1)];
line_ally=[line_top(:,2);line_right(:,2);line_bottom(:,2);line_left(:,2)];
line_allz=line_ally*0+1;

%rotate these points around the sphere
line_r = [line_allx(:), line_ally(:), line_allz(:)].*1./sqrt((line_allx(:).^2+line_ally(:).^2+1));
line_r2 = line_r*rotmat';

%convert these points to the stereogram
line_r_lambda=1./(line_r2(:,3)+1);

%convert from r_exp to the stereogram
line_P_rx=line_r2(:,1).*line_r_lambda;
line_P_ry=line_r2(:,2).*line_r_lambda;

%plot the blue box
plot(line_P_rx,line_P_ry,'b','LineWidth',2);

%
%do the fundamental zone
%write three axes, and construct the lines that go between one and the next
ax1=[0 0 1];
ax2=[1 1 2];
ax3=[1 0 1];

%normalize length
ax1=ax1./norm(ax1);
ax2=ax2./norm(ax2);
ax3=ax3./norm(ax3);
nl=100;

%inscribe along the great circles and equispace
ang_12=acos(dot(ax1,ax2));
ax_12=cross(ax1,ax2);
ax_12=ax_12/norm(ax_12);
ax_2p=cross(ax_12,ax1);
ax_2p=ax_2p/norm(ax_2p);
sample_spacing12=linspace(0,ang_12,nl);
line_12=ax1'*cos(sample_spacing12)+ax_2p'*sin(sample_spacing12);

ang_23=acos(dot(ax2,ax3));
ax_23=cross(ax2,ax3);
ax_23=ax_23/norm(ax_23);
ax_3p=cross(ax_23,ax2);
ax_3p=ax_3p/norm(ax_3p);
sample_spacing23=linspace(0,ang_23,nl);
line_23=ax2'*cos(sample_spacing23)+ax_3p'*sin(sample_spacing23);

ang_31=acos(dot(ax3,ax1));
ax_31=cross(ax3,ax1);
ax_31=ax_31/norm(ax_31);
ax_1p=cross(ax_31,ax3);
ax_1p=ax_1p/norm(ax_1p);
sample_spacing31=linspace(0,ang_31,nl);
line_31=ax3'*cos(sample_spacing31)+ax_1p'*sin(sample_spacing31);

%projected on stereogram
line_set=[line_12';line_23';line_31'];
line_set_lambda=1./(line_set(:,3)+1);
line_set_rx=line_set(:,1).*line_set_lambda;
line_set_ry=line_set(:,2).*line_set_lambda;

plot(line_set_rx(1:end),line_set_ry(1:end),'y', 'LineWidth',2);


% %plot the individual segments
% plot(line_set_rx(1:end),line_set_ry(1:end),'r');
% scatter(line_set_rx(201),line_set_ry(201),'m');
% plot(line_set_rx(nl+(1:nl)),line_set_ry(nl+(1:nl)),'g');
% plot(line_set_rx(2*nl+(1:nl)),line_set_ry(2*nl+(1:nl)),'b');

%% lets inflate this up to the full pattern
nline=20;
line_P_rxn=zeros(4*nline,24);
line_P_ryn=zeros(4*nline,24);
stereo_remap_maskp=zeros(num_stereo,num_stereo,24);
stereo_remap_2dnp=stereo_remap_maskp;

cs_al=loadCIF('Al-Aluminum.cif');

EBSP_av=EBSD_simulation;
[Pat_sim_astroind,r_exp]= EBSP_gen( EBSD_simulation,rotmat,screen_int); %generate the EBSP for this iteration

r_n=zeros(size(r_exp,1),size(r_exp,2),24);
n_sym=24;
for n=1:n_sym
    sym_matrix=cs_al.rot(n).matrix;
    r_n(:,:,n) = r_exp*sym_matrix;
end

stereo_x_ok=stereo_x(:);
stereo_y_ok=stereo_y(:);

% calculate and plot each symmetrically equivalent pattern in the
% stereogram
parfor p=1:24

    x_lookup=EBSP_av.xpts_screen;
    y_lookup=EBSP_av.ypts_screen;

    x_lookup=x_lookup(:,2:end-1);
    y_lookup=y_lookup(:,2:end-1);
    pat1_plot_n=pat1_plot(:,2:end-1);

    x_lookup=x_lookup(2:end-1,1);
    y_lookup=y_lookup(2:end-1,1);
    pat1_plot_n=pat1_plot_n(2:end-1,:);

    r = [x_lookup(:), y_lookup(:), y_lookup(:)*0+1].*1./sqrt((x_lookup(:).^2+y_lookup(:).^2+1));
    r_exp2 = r*rotmat';

    r_3=r_n(:,3,p);
    r_3(r_3<0)=-r_3(r_3<0);
    r_lambda=1./(r_3+1);

    %convert from r_exp to the stereogram
    P_rxn=r_n(:,1,p).*r_lambda;
    P_ryn=r_n(:,2,p).*r_lambda;

    % interpolate
    stero_interp_exp_n=scatteredInterpolant(P_rxn(:),P_ryn(:),pat1_plot(:),'natural','none');

    %plot the outline of the detector
    line_top=[EBSD_simulation.xpts_screen(:,1),EBSD_simulation.ypts_screen(:,1)];
    line_bottom=flipud([EBSD_simulation.xpts_screen(:,end),EBSD_simulation.ypts_screen(:,end)]);
    line_right=[EBSD_simulation.xpts_screen(end,:);EBSD_simulation.ypts_screen(end,:)]';
    line_left=flipud([EBSD_simulation.xpts_screen(1,:);EBSD_simulation.ypts_screen(1,:)]');
    
    s2=linspace(1,254,nline);
    s2(2:end-1)=round(s2(2:end-1));

    line_allx=[line_top(s2,1);line_right(s2,1);line_bottom(s2,1);line_left(s2,1)];
    line_ally=[line_top(s2,2);line_right(s2,2);line_bottom(s2,2);line_left(s2,2)];
    line_allz=line_ally*0+1;

    %rotate these points around the sphere
    line_r = [line_allx(:), line_ally(:), line_allz(:)].*1./sqrt((line_allx(:).^2+line_ally(:).^2+1));
    line_r2 = line_r*rotmat'*cs_al.rot(p).matrix;
    line_r2_3=line_r2(:,3);
    line_r2_3(line_r2_3<0)=-line_r2_3(line_r2_3<0);
    %convert these points to the stereogram
    line_r_lambda=1./(line_r2_3+1);

    %convert from r_exp to the stereogram
    line_P_rx=line_r2(:,1).*line_r_lambda;
    line_P_ry=line_r2(:,2).*line_r_lambda;
    
    %create the blank frame
    stereo_remap_2dn=zeros(num_stereo,num_stereo);
    stereo_remap_mask=ones(num_stereo,num_stereo);
num_el=numel(stereo_remap_mask);
    % %find the points in the detector frame
    % [in,on] = inpolygon(stereo_x_ok,stereo_y_ok,line_P_rx,line_P_ry);
    in=true(num_el,1);
    stereo_remap_mask_0=stereo_remap_mask*0;
    %read the intepolant
    stereo_remap_n=stero_interp_exp_n(stereo_x_ok(in),stereo_y_ok(in));
    stereo_remap_2dn(in)=stereo_remap_n;
    stereo_remap_mask(in)=stereo_remap_mask_0(in);

    % stereo_remap_2dn=reshape(stereo_remap_n,num_stereo,num_stereo);
    % stereo_mask_2dn=isnan(stereo_remap_2dn);
    
    % Figure 6,7
    %plot the base dynamical pattern
    figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
    axis image; axis xy; colormap('gray');
    hold on;

    %plot the remapped experiment
    i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2dn);
    axis image; axis xy; colormap('gray');
    i1.AlphaData=1-stereo_remap_mask;

    %plot the blue box
    plot(line_P_rx,line_P_ry,'b');
    xlim([-1 1]);
    ylim([-1 1]);

    title(int2str(p));

    line_P_rxn(:,p)=line_P_rx;
    line_P_ryn(:,p)=line_P_ry;
    stereo_remap_maskp(:,:,p)=stereo_remap_mask;
    stereo_remap_2dnp(:,:,p)=stereo_remap_2dn;
end

%% Individual plots of the 24 equivalent patterns reprojected 
% (this is commented out by default to reduce window count)

% for p=1:24
%     %plot the base dynamical pattern
%     figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
%     axis image; axis xy; colormap('gray');
%     hold on;
% 
%     %plot the remapped experiment
%     i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2dnp(:,:,p));
%     axis image; axis xy; colormap('gray');
%     i1.AlphaData=1-stereo_remap_maskp(:,:,p);
% 
%     %plot the blue box
%     plot(line_P_rxn(:,p),line_P_ryn(:,p),'b');
%     xlim([-1 1]);
%     ylim([-1 1]);
%     title(int2str(p));
% end

%% Plot all outlines of the reprojected experimental patterns 
% on a simulated pattern

figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray');
hold on;

for p=1:24
    %plot the blue box
    plot(line_P_rxn(:,p),line_P_ryn(:,p),'b', 'LineWidth',1.5);
    xlim([-1 1]);
    ylim([-1 1]);
    title(int2str(p));
end

%% Plot the sterogram and try to work out the points for this projection

%read the experimental pattern
pat1_plot2=pat1_plot;
pat_c=pat1_plot2(:);

pat_maps=zeros(2001,2001,24);
pat_mapn=pat_maps;

%pick the symemtry
for p=1:24
    %extract the points of this pattern
    %rotate these points of the detector around the sphere
    box_x=EBSD_simulation.xpts_screen(:);
    box_y=EBSD_simulation.ypts_screen(:);
    box_z=box_y*0+1;

    % box_r=[box_x,box_y,box_z];
    box_r = [box_x(:), box_y(:), box_z(:)].*1./sqrt((box_x(:).^2+box_y(:).^2+1));

    box_r2=box_r*rotmat'*cs_al.rot(p).matrix;
    box_r2x=box_r2(:,1);
    box_r2y=box_r2(:,2);
    box_r2z=box_r2(:,3);

    %now find those which are in the positive and negative hemisphere ('p'
    %and 'n' points) based upon their z coords
    box_r2z_pz=find(box_r2z>0);
    box_r2z_nz=find(box_r2z<0);

    %subselect the XYZ positives
    box_r2xp=box_r2x(box_r2z_pz);
    box_r2yp=box_r2y(box_r2z_pz);
    box_r2zp=box_r2z(box_r2z_pz);

    %subselect the XYZ negatives
    box_r2xn=-box_r2x(box_r2z_nz);
    box_r2yn=-box_r2y(box_r2z_nz);
    box_r2zn=-box_r2z(box_r2z_nz);

    %extract the pattern data (we are going to use scatter for a lazy plot)
    pat_cp=pat_c(box_r2z_pz);
    pat_cn=pat_c(box_r2z_nz);

    %sterographic projection
    box_r2_lambdap=1./(box_r2zp+1);
    box_r2_xp=box_r2xp.*box_r2_lambdap;
    box_r2_yp=box_r2yp.*box_r2_lambdap;

    %find the boundary of the point cloud
    boundary_r2_p=boundary(box_r2_xp,box_r2_yp);
    boundary_r2_px=box_r2_xp(boundary_r2_p);
    boundary_r2_py=box_r2_yp(boundary_r2_p);

    %sterographic projection
    box_r2_lambdan=1./(box_r2zn+1);
    box_r2_xn=box_r2xn.*box_r2_lambdan;
    box_r2_yn=box_r2yn.*box_r2_lambdan;

    %find the boundary of the point cloud
    boundary_r2_n=boundary(box_r2_xn,box_r2_yn);
    boundary_r2_nx=box_r2_xn(boundary_r2_n);
    boundary_r2_ny=box_r2_yn(boundary_r2_n);
  
    %find the pounts within each boundary
    sx=stereo_x;
    sy=stereo_y;

    [in_p]=inpolygon(sx,sy,boundary_r2_px(:),boundary_r2_py(:));
    [in_n]=inpolygon(sx,sy,boundary_r2_nx(:),boundary_r2_ny(:));

    x_lookup=EBSD_simulation.xpts_screen;
    y_lookup=EBSD_simulation.ypts_screen;

    r = [x_lookup(:), y_lookup(:), y_lookup(:)*0+1].*1./sqrt((x_lookup(:).^2+y_lookup(:).^2+1));
    r_exp2 = r*rotmat'*cs_al.rot(p).matrix;
    %create the interpolant
    r_3=r_exp2(:,3);
    r_lambda=1./(r_3+1);

    %convert from r_exp to the stereogram
    P_rxn=r_exp2(:,1).*r_lambda;
    P_ryn=r_exp2(:,2).*r_lambda;

    % generate interpolate
    stero_interp_exp_s=scatteredInterpolant(P_rxn(:),P_ryn(:),pat1_plot(:),'natural','none');

    stero_interp_exp_sp=scatteredInterpolant(box_r2_xp(:),box_r2_yp(:),pat_cp(:),'natural','none');
    stero_interp_exp_sn=scatteredInterpolant(box_r2_xn(:),box_r2_yn(:),pat_cn(:),'natural','none');

    P_stero_xp=sx(in_p);
    P_stero_yp=sy(in_p);
    Pat_interp_p=stero_interp_exp_sp(P_stero_xp,P_stero_yp);

    P_stero_xn=sx(in_n);
    P_stero_yn=sy(in_n);
    Pat_interp_n=stero_interp_exp_sn(P_stero_xn,P_stero_yn);

    pat_map=0*sx;
    pat_num=0*sx;
    pat_num(in_p)=pat_num(in_p)+1;
    pat_num(in_n)=pat_num(in_n)+1;

    pat_map(in_p)=Pat_interp_p;
    pat_map(in_n)=Pat_interp_n;

    pat_maps(:,:,p)=pat_map;
    pat_mapn(:,:,p)=pat_num;
end

%% Normalize the overlaid stereogram and plot
pat_maps_num=sum(pat_mapn,3);
pat_maps_sum=sum(pat_maps,3);
pat_maps_norm=pat_maps_sum./pat_maps_num;

%% final stereogram
figure; 
subplot(1,3,1);
imagesc(stereo_xl,stereo_yl,pat_maps_sum); % total stereogram
axis image; axis xy; colormap('gray'); axis off;
subplot(1,3,2);
imagesc(stereo_xl,stereo_yl,pat_maps_norm); % stereogram normalized by weight
axis image; axis xy; colormap('gray'); axis off;
% sum/weight hemisphere
subplot(1,3,3); imagesc(stereo_xl,stereo_yl,pat_maps_num); % overlapped patterns
axis image; axis xy; colormap('gray'); axis off; 

%% Difference between simulated and reprojected patterns
% normalize the reprojected dynamics pattern
pat_maps_dyn=reshape(stereo_pat_plot,num_stereo,num_stereo);
pat_maps_dynn=(pat_maps_dyn-mean(pat_maps_dyn(:)))./std(pat_maps_dyn(:));
pat_maps_normn=(pat_maps_norm-mean(pat_maps_norm(:)))./std(pat_maps_norm(:));

figure; imagesc(stereo_xl,stereo_yl,pat_maps_norm-pat_maps_dyn);
axis image; axis xy; colormap('gray');
hold on;
clim([-6 6]);

%% Write some sterograms to tifs


imwrite(uint16(flipud(normalizeto16bit(pat_maps_norm))), fullfile(path1,'results','stereo_stack_exp.tif'));
imwrite(uint16(flipud(normalizeto16bit(pat_maps_dyn))), fullfile(path1,'results','stereo_stack_dyn.tif'));
     

%%
line_r = [line_allx(:), line_ally(:), line_allz(:)].*1./sqrt((line_allx(:).^2+line_ally(:).^2+1));
line_r2 = line_r*rotmat'*cs_al.rot(p).matrix;
line_r2_3=line_r2(:,3);

line_p=line_r2_3*0;
line_y1=line_r2_3*0+1;
line_p(line_r2_3<0)=line_y1(line_r2_3<0);
line_p=logical(line_p);
line_n=logical(1-line_p);

line_r2_3(line_r2_3<0)=-line_r2_3(line_r2_3<0);


%convert these points to the stereogram
line_r_lambda=1./(line_r2_3+1);

%convert from r_exp to the stereogram
line_P_rx=line_r2(:,1).*line_r_lambda;
line_P_ry=line_r2(:,2).*line_r_lambda;

%plot the base dynamical pattern
figure; imagesc(stereo_xl,stereo_yl,reshape(stereo_pat_plot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray');
hold on;

%plot the remapped experiment
i1=imagesc(stereo_xl,stereo_yl,stereo_remap_2dnp(:,:,p));
axis image; axis xy; colormap('gray');
i1.AlphaData=1-stereo_remap_maskp(:,:,p);

%plot the blue box
plot(line_P_rx,line_P_ry,'b');

%plot the positive
line_P_rx_p=line_P_rx(line_p);
line_P_ry_p=line_P_ry(line_p);
plot(line_P_rx_p,line_P_ry_p,'y');

%plot the negative
line_P_rx_n=line_P_rx(line_n);
line_P_ry_n=line_P_ry(line_n);
plot(line_P_rx_n,line_P_ry_n,'m');

[in_l]=inpolygon(line_P_rx_n,line_P_ry_n,line_P_rx_p,line_P_ry_p);
line_P_rx_n1=line_P_rx_n(in_l);
line_P_ry_n1=line_P_ry_n(in_l);

line_P_rx_n1_r=sqrt(line_P_rx_n1.^2+line_P_ry_n1.^2);

line_P_rx_n1=line_P_rx_n1./line_P_rx_n1_r;
line_P_ry_n1=line_P_ry_n1./line_P_rx_n1_r;
plot(line_P_rx_n1,line_P_ry_n1,'w','linewidth',3);

xlim([-1 1]);
ylim([-1 1]);
title(int2str(p));

% plot a circle in green
thetac=linspace(0,2*pi,360);
x_c=1*sin(thetac);
y_c=1*cos(thetac);
plot(x_c,y_c,'g')