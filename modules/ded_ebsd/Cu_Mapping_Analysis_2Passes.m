%% Cu_Mapping_Analysis_2Passes.m
% By T Ben Britton, November 2023
% Edited and annotated by Tianbi Zhang, November 2023

% This is the script for the 20x20 Cu mapping dataset.
% The first pass performs a PC search for every pattern.
% PC of the 1st pattern of the series (theoretically any one) should
% be known for the PC fit model in the first pass to work.
% PCx, PCy and PCz values are then fitted to sep[arate linear functions of
% beam shift in X and Y to build and correct the previous PC model based on 
% assumptions of geometries.
% The second indexing pass is then performed using the new PC data.

% Requirements:
% Cu patterns (see doi:10.5281/zenodo.8342489)
% MTEX toolbox
% MATLAB image processing, curve fitting and statisics and machine learning
% toolboxes

%% Start of the script
close all;
clear;
home;

%% start up Astro etc.
InputUser.location_mtex="C:\Users\billy\Documents\MATLAB\mtex-5.10.0";
run(fullfile(InputUser.location_mtex, "startup_mtex.m"));
InputUser.Astro_loc = "C:\Users\billy\Documents\MATLAB\AstroEBSD-main";
run(fullfile(InputUser.Astro_loc, "start_AstroEBSD.m"));

%% location of patterns
pattern_loc='C:\Users\billy\OneDrive - UBC\PhD\TKD\20230807_Cu_EBSD';
pattern_root='Spot';

%% AstroEBSD Parameters
% BG correction
Settings_Cor.radius=1;
Settings_Cor.radius_frac=0.85;
Settings_Cor.gfilt=1; %use a low pass filter
Settings_Cor.gfilt_s = 10; %low pass filter sigma - you can change this

%single pattern
Settings_PCin.start=[0.52,0.45,0.605];
Settings_PCin.range=[0.05 0.05 0.05]; %+- these values

%Set up the radon transform peak finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak hunt
Settings_Rad.max_peaks=10; %max number of peaks to return
Settings_Rad.num_peak=35; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=10; %search size in theta steps
Settings_Rad.rho_search_per=0.15; %radon search in fractions
Settings_Rad.min_peak_width=0.004; %min rseperation of the peak width, in pixels

%build the phases
InputUser.Phase_Folder = fullfile(InputUser.Astro_loc,'phases\phasefiles');
InputUser.Phase_Input  = {'Cu'}; %Si, Ferrite
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

% For the first pass, we assume the sample tilt is 20 degrees.
Detector_tilt = RTM.Rx(-20 * pi / 180);

%% Load the EBSPs
disp('Loading EBSPs into RAM')
pattern_all=zeros(252,252,400); 
for pat=1:400
    pattern_name=[pattern_root int2str(pat)];
    pattern_read = hdfpatternread(pattern_loc, pattern_name, true);
    [pattern_all(:,:,pat),~] = EBSP_BGCor( pattern_read,Settings_Cor);
end

%% First indexing pass (with PC refinement)

tic
for pat=1:400
    pat

EBSP_One=struct;
pattern_read=pattern_all(:,:,pat);
[EBSP_One.PatternIn,Settings_Cor ] = EBSP_BGCor( pattern_read,[]);

%find the PC
[EBSP_One_PC_S] = EBSP_PCSearch(EBSP_One,Settings_Cor,Settings_Rad,Settings_PCin,Phase_Num,Crystal_LUT,Crystal_UCell);
%Index normally
[ EBSP_One.Peak_Centre,EBSP_One.Single.Peak_Set_All,EBSP_One.Peak_Set_All,...
            EBSP_One.R_EBSP,EBSP_One.R_Edge,EBSP_One.R_rho,EBSP_One.R_theta ] ...
            = EBSP_RadHunt( EBSP_One.PatternIn,Settings_Rad);

% Convert the bands to normal space
[ EBSP_One.nhat_gnom] = EBSP_NormConv( EBSP_One.Peak_Centre,size(EBSP_One.PatternIn),EBSP_One_PC_S.PC_out);

% Index
[EBSP_One.rotdata{1},EBSP_One.banddata]=EBSP_Index(EBSP_One.nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},Detector_tilt); %#ok<PFBNS>

%generate the geometry
% [ EBSP_One.PatternGeometry ] = EBSP_Gnom( Settings_Cor,EBSP_One_PC_S.PC_out );
toc
% EBSP_OneFigure=Plot_SinglePattern(EBSP_One,Crystal_UCell,Crystal_LUT,1);

s_err_old(pat)=EBSP_One.rotdata{1}.error;
EBSP_One_all{pat}.Pattemetry.PC=EBSP_One_PC_S.PC_out;
end

%% check an index - uncomment if you need/want
% n=243;
% EBSP_One=EBSP_One_all{n};
% Astro_EBSPset(EBSP_One,Settings_Cor,Settings_Rad,Settings_PCin,InputUser);

%% extract the PC values
for n=1:400; PC(:,n)=EBSP_One_all{n}.Pattemetry.PC; end

PC_x=reshape(PC(1,:),20,20)';
PC_y=reshape(PC(2,:),20,20)';
PC_z=reshape(PC(3,:),20,20)';

figure; 
subplot(1,3,1);
imagesc(PC_x);
axis image; axis xy; axis tight; colormap('gray');
subplot(1,3,2);
imagesc(PC_y);
axis image; axis xy; axis tight; colormap('gray');
subplot(1,3,3);
imagesc(PC_z);
axis image; axis xy; axis tight; colormap('gray');


%% Fit the PC data
pos_vals=1:20;
pos_vals=pos_vals-mean(pos_vals);
[xmap,ymap]=meshgrid(pos_vals,pos_vals);

% Set up fittype and options.
ft = fittype( 'poly11' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';

[xData_PCx, yData_PCx, zData_PCx] = prepareSurfaceData( xmap, ymap, PC_x );
[xData_PCy, yData_PCy, zData_PCy] = prepareSurfaceData( xmap, ymap, PC_y );
[xData_PCz, yData_PCz, zData_PCz] = prepareSurfaceData( xmap, ymap, PC_z );

[fitresult_PCx, gof_PCx] = fit( [xData_PCx, yData_PCx], zData_PCx, ft, opts );
[fitresult_PCy, gof_PCy] = fit( [xData_PCy, yData_PCy], zData_PCy, ft, opts );
[fitresult_PCz, gof_PCz] = fit( [xData_PCz, yData_PCz], zData_PCz, ft, opts );

PCx=fitresult_PCx.p00;
PCy=fitresult_PCy.p00;
PCz=fitresult_PCz.p00;

PC_model=[PCx,PCy,PCz]; %central value from the model

%extract the map data
sample_tilt=atand(fitresult_PCy.p01/fitresult_PCz.p01);
step_size=fitresult_PCx.p10*55*251; %55 um per pixel, and 251 pixels per detector

%{
Scan type = regular scan over a spaced grid
X start: 600, X spacing=10, Y start: 100, Y spacing=10
Grid points: 20x20
Mag: 112.0083x
Pixel size (m): 1.005117E-06
Stage X (m): 0.09806616
Stage Y (m): 0.06753983
Stage Z (m): 0.0009811198
%}

step_sem=1.005117E-06*10*1E6; %spacing * pixel size * 1E6 to give in um

%% Now we can generate a model of the PC from this data and plot this

PC_x_map=PC_model(1)+xmap*fitresult_PCx.p10;
PC_y_map=PC_model(2)+ymap*fitresult_PCy.p01;
PC_z_map=PC_model(3)+ymap*fitresult_PCz.p01;

figure; 
subplot(2,3,1);
imagesc(PC_x);
axis image; axis xy; axis tight; colormap('gray'); clim([min(PC_x_map(:)) max(PC_x_map(:))]);
subplot(2,3,2);
imagesc(PC_y);
axis image; axis xy; axis tight; colormap('gray'); clim([min(PC_y_map(:)) max(PC_y_map(:))]);
subplot(2,3,3);
imagesc(PC_z);
axis image; axis xy; axis tight; colormap('gray'); clim([min(PC_z_map(:)) max(PC_z_map(:))]);
subplot(2,3,4);
imagesc(PC_x_map);
axis image; axis xy; axis tight; colormap('gray');
subplot(2,3,5);
imagesc(PC_y_map);
axis image; axis xy; axis tight; colormap('gray');
subplot(2,3,6);
imagesc(PC_z_map);
axis image; axis xy; axis tight; colormap('gray');

%% Second indexing pass using the new PC model
%update max peaks to use to make this indexing more robust
Settings_Rad.max_peaks=7;

pat_nums=1:400;
pat_nums_square=reshape(pat_nums,20,20)';

Settings_Cor.radius=1;
Settings_Cor.radius_frac=0.85;

%Define all rotation matrices needed in the code
Detector_tilt = RTM.Rx(-(90-sample_tilt) * pi / 180);

tic
for n=1:400
    EBSP_One=struct;
    pattern_read=pattern_all(:,:,pat_nums_square(n));
    [EBSP_One.PatternIn,Settings_Cor ] = EBSP_BGCor( pattern_read,Settings_Cor);
    EBSP_One.PC=[PC_x_map(n),PC_y_map(n),PC_z_map(n)];
    
    %Index normally
    [ EBSP_One.Peak_Centre,EBSP_One.Single.Peak_Set_All,EBSP_One.Peak_Set_All,...
        EBSP_One.R_EBSP,EBSP_One.R_Edge,EBSP_One.R_rho,EBSP_One.R_theta ] ...
        = EBSP_RadHunt( EBSP_One.PatternIn,Settings_Rad);

    % Convert the bands to normal space
    [ EBSP_One.nhat_gnom] = EBSP_NormConv( EBSP_One.Peak_Centre,size(EBSP_One.PatternIn),EBSP_One.PC);

    % Index
    [EBSP_One.rotdata{1},EBSP_One.banddata]=EBSP_Index(EBSP_One.nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},Detector_tilt); %#ok<PFBNS>
    EBSP_One_re{n}=EBSP_One;
    n
    toc
end

%% record and check the error values and Euler Angles. 
% you can do the same for the first pass too.
for n=1:400
X(n) = mod(n-1,20) + 1;
Y(n) = floor((n-1)/20) + 1;
EulerZ1(n) = EBSP_One_re{n}.rotdata{1}.eang(1);
EulerX(n) = EBSP_One_re{n}.rotdata{1}.eang(2);
EulerZ2(n) = EBSP_One_re{n}.rotdata{1}.eang(3);
s_err(n)=EBSP_One_re{n}.rotdata{1}.error;
end
figure;
subplot(1,2,1);plot(s_err);

subplot(1,2,2);
imagesc(reshape(s_err*180/pi,20,20)');
caxis([0 1E-2]*180/pi); axis image; axis xy;

%% check a pattern index (optional)
% pat_check=10;
% EBSP_One=EBSP_One_re{pat_check};
% [ EBSP_One.PatternGeometry ] = EBSP_Gnom( Settings_Cor,EBSP_One.PC);
% 
% EBSP_OneFigure=Plot_SinglePattern(EBSP_One,Crystal_UCell,Crystal_LUT,1);
% Settings_PCin.start=EBSP_One.PC;
% Astro_EBSPset(EBSP_One,Settings_Cor,Settings_Rad,Settings_PCin,InputUser);

%% Create stuff in MTEX
cs = loadCIF('C:\Users\billy\Documents\GitHub\AstroEBSD_v2\phases\cifs\Cu-copper.cif');
% Pat 1 is from top left because X sense of the beam is the opposite to X
% of the detector.
setMTEXpref('xAxisDirection','east'); 
setMTEXpref('zAxisDirection','intoplane');

% build the coordinate maps. To get the accurate micron bar, we multiply by
% the step size.
prop.y = double(X')*10; 
prop.x = double(Y')*10;
prop.MAE=double(s_err);

ori = ...
  rotation('Euler',EulerZ1,EulerX,EulerZ2);

try
ebsd = EBSD(ori.', ones(size(ori)),{'notIndexed',cs},'Options',prop);
catch
    ebsd = EBSD(ori.', ones(size(ori)),{'notIndexed',cs},prop);
end

ipfKey = ipfHSVKey(ebsd('copper'));
figure;
plot(ipfKey);

colorKey = ipfHSVKey(cs);
colorKey.inversePoleFigureDirection = xvector;
color = colorKey.orientation2color(ebsd('indexed').orientations);

cShape = crystalShape.cube(ebsd.CS);

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'));
big_grains = grains(grains.grainSize > 10);

figure;
plot(ebsd,colorKey.orientation2color(ebsd.orientations));
hold on;
plot(big_grains('copper'),0.4*cShape,'linewidth',2,'colored','micronbar','off');
hold off;
legend off;

% You can compare the errors using the old and new PC methods.
figure;
histogram(s_err_old,'EdgeColor','none','Normalization','count');
hold on;
histogram(s_err,'EdgeColor','none','Normalization','count');
xlabel("Mean Angular Error (\circ)");
ylabel("Counts");
