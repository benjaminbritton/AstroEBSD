clear
home 
close all

% This script was created to kick off discussions at the 2019 EBSD meeting
% Deck to show how pattern matching can work with MTEX

%% Toolbox locations for AstroEBSD and MTEX

%these are Ben's settings
InputUser.Astro_loc='C:\Users\tbritton\Documents\GitHub\AstroEBSD_v2\'; %Change this to your AstroEBSD location
%  InputUser.Astro_loc='C:\Users\bbrit\Documents\GitHub\AstroEBSD_v2'; %Change this to your AstroEBSD location

%  location_mtex='C:\Users\bbrit\Documents\mtex-5.4.0'; %Change this to where you keep your MTEX folder
location_mtex='E:\Communal_MatlabPlugins\mtex-5.4.0';

run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m'));
run(fullfile(location_mtex,'startup_mtex.m'));


%% Load the experimental data we want to work with
InputUser.Phase_Folder = fullfile(InputUser.Astro_loc,'phases\phasefiles'); %location of the pha files
% InputUser.HDF5_folder='C:\Users\bbrit\Documents\Issues\BigFile\'; %Change this to the file location in whch you have saved the example data
InputUser.HDF5_folder='E:\Ben\';

InputUser.HDF5_file= 'Demo_Ben.h5';
InputUser.Phase_Input  = {'Ferrite'};


%% Set up the RTM

%setttings for RTM
RTM.screensize = 128; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM.iterations = 4;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 128; %LPT size used in pixels

%From AstroEBSD
%background correction
Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=5; %low pass filter sigma
Settings_CorX.radius=0; %use a radius mask
Settings_CorX.radius_frac=0.85; %fraction of the pattern width to use as the mask
Settings_CorX.hotpixel=0; %hot pixel correction
Settings_CorX.hot_thresh=1000; %hot pixel threshold
Settings_CorX.resize=1; %resize correction
Settings_CorX.size=RTM.screensize; %image height
Settings_CorX.RealBG=0; %use a real BG
Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG
Settings_CorX.SquareCrop = 1; %make square the EBSP
Settings_CorX.SplitBG=1; %deal with a split screen

%% Low level setting up stuff - you shouldn't need to change this
RTM.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'masterpatterns'); %location of the binary files used for RTM

[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation


%% Read data
[ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadHDF5( InputUser );

%load the data from h5
[MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadHDF5( InputUser );
[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);

%adjust the pattern centre for the square cropping - this corrects PCx
[ MapInfo.Data_InputMap ] = PC_square( MapInfo.EBSPData, MapInfo.Data_InputMap,Settings_CorX );

%get stuff ready for XCF
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

%% create a reference pattern from the Dynamics bin file
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},RTM.Phase_Folder);
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%% Plot the sphere (because we can)
%create a spherical pattern
[rx,ry,rz]=sphere(1000);
%if you rotate rx, ry and rz you can rotate the pattern
[i_data] = Cube_Sample(rx(:),ry(:),rz(:),screen_int,RTM_info.isHex);

figure;
surf(rx,ry,rz,reshape(i_data,size(rx,1),size(rx,2)),'EdgeColor','none'); 
colormap('gray');
axis off; axis equal;
axis vis3d;

%% Generate a simulated pattern
% create a geometry
PC_simulation=[0.5,0.5,0.4]; %try varying this
[EBSD_simulation ] = EBSP_Gnom( RTM,PC_simulation); %you can change PC_in if you want

%create an orientation
eangs=[25 20 0]*pi/180; %phi1, PHI, phi2
gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
[ Pat_sim_eang ] = EBSP_gen( EBSD_simulation,gmatrix*Detector_tilt,screen_int); %generate the EBSP for this iteration 

figure;
subplot(1,1,1);
I1=pPattern(Pat_sim_eang,EBSD_simulation);
title('Simulated pattern');

%% Now generate a pattern from the Bruker orientation data

P_test=10; %the 10th pattern captured
[P_r,P_c]=fPointFind(P_test,MapInfo.Data_InputMap.PMap);

%read the experimental pattern
[ Pat_exp ] = bReadEBSP(EBSPData,P_test);
[ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,Settings_CorX);

%Read the PC from Bruker data
PC_in=[MapInfo.Data_InputMap.PCX(P_r,P_c),MapInfo.Data_InputMap.PCY(P_r,P_c),MapInfo.Data_InputMap.DD(P_r,P_c)];
[EBSD_geom ] = EBSP_Gnom( RTM,PC_in); %you can change PC_in if you want

%read the Euler angles
eangs=[MapInfo.Data_InputMap.phi1(P_r,P_c),MapInfo.Data_InputMap.PHI(P_r,P_c),MapInfo.Data_InputMap.phi2(P_r,P_c)]*pi/180;
gmatrix=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));

%simulate
[ Pat_sim_B] = EBSP_gen( EBSD_geom,gmatrix*Detector_tilt,screen_int); 

figure;
subplot(1,2,1);
I1=pPattern(Pat_exp_cor,EBSD_geom);
title('Experimental pattern');

subplot(1,2,2);
I1=pPattern(Pat_sim_B,EBSD_geom);
title('Simulated pattern');

%% Refine the Bruker solution

%prepare the experimental pattern for refinement
[Pat_Ref_r,XCF_data_fill] = refine_prep(Pat_exp_cor,SettingsXCF,RTM);

%refine
[G_Refined,regout_R] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,gmatrix*Detector_tilt,SettingsXCF,screen_int,RTM_info.isHex,RTM);
%[orientation,XCF values]=refine5(Pattern,mean-geometry,updatedPC,orientation,SettingsXCF,screen_int,isHex,RTM)
%XCF values = [x shift, y shift, ~ , normalised chi]

%update simulation
[ Pat_sim_B_ref ] = EBSP_gen( EBSD_geom,G_Refined,screen_int,screen_int.isHex ); %generate the EBSP for this iteration 

figure;
subplot(2,2,1);
I1=pPattern(Pat_exp_cor,EBSD_geom);
title('Experimental pattern');

subplot(2,2,2);
I1=pPattern(Pat_sim_B-Pat_sim_B_ref,EBSD_geom);
title('Difference');

subplot(2,2,3);
I1=pPattern(Pat_sim_B,EBSD_geom);
title('Refined Simulation pattern');

subplot(2,2,4);
I1=pPattern(Pat_sim_B_ref,EBSD_geom);
title('Bruker Simulation pattern');

%% Now do SO(3) searching
%populate fundamental zone of SO(3)
cs_phase=loadCIF(RTM_info.cif_file);
[ library_G ] = SO3_rotmat_gen( cs_phase,RTM.Sampling_Freq);

%generate the library
[ template_library ] = Library_Gen(EBSD_geom,screen_int,RTM_info.isHex,library_G,2,SettingsXCF);

%compare with the library
[G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r.FFT,SettingsXCF,XCF_data_fill,1);

%refine
[G_Refined_SO3,regout_R] = refine5(Pat_Ref_r,EBSD_geom,EBSD_geom.PC,G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);

%update simulation
[ Pat_sim_SO3 ] = EBSP_gen( EBSD_geom,G_SO3,screen_int,screen_int.isHex ); %generate the EBSP for this iteration 
%update simulation
[ Pat_sim_SO3r ] = EBSP_gen( EBSD_geom,G_Refined_SO3,screen_int,screen_int.isHex ); %generate the EBSP for this iteration 

figure;
subplot(2,2,1);
I1=pPattern(Pat_exp_cor,EBSD_geom);
title('Experimental pattern');

subplot(2,2,2);
I1=pPattern(Pat_sim_SO3-Pat_sim_SO3r,EBSD_geom);
title('Difference');

subplot(2,2,3);
I1=pPattern(Pat_sim_SO3,EBSD_geom);
title('SO(3) Template Found');

subplot(2,2,4);
I1=pPattern(Pat_sim_SO3r,EBSD_geom);
title('SO(3) Refined');

%% Prepare data from a h5 to run as map

pat_list=1:MicroscopeData.NPoints;
%read the the patterns and prepare them for XCF

%prepare the data fill for the faster FFT
[~,XCF_data_fill] = refine_prep(zeros(RTM.screensize),SettingsXCF,RTM);

%create empty arrays in the structure
Pat_Ref_r = struct('FFT', zeros(numel(XCF_data_fill), numel(XCF_data_fill)), 'logp', zeros(RTM.LPTsize,RTM.LPTsize),'pat_in',zeros(RTM.screensize));

%build the structure with the input data
parfor p=1:numel(pat_list)
    [ Pat_exp ] = bReadEBSP(EBSPData,p);
    [ Pat_exp_cor ] = EBSP_BGCor( Pat_exp,Settings_CorX);
    [ Pat_Ref_r(p) ] = refine_prep(Pat_exp_cor,SettingsXCF,RTM); %note this stores every pattern, the FFT and the log-polar 
end

%% Extract useful mapdata
map_pcx=MapInfo.MapData.PCX;
map_pcy=MapInfo.MapData.PCY;
map_pcz=MapInfo.MapData.DD;

map_phi1=map_pcx;
map_PHI=map_pcx;
map_phi2=map_pcx;

g_bruker=zeros(3,3,numel(map_pcx));
g_bruker_det=zeros(3,3,numel(map_pcx));

for p=1:numel(map_pcx)
    %find the point in the map
    [P_r,P_c]=fPointFind(p,MapInfo.Data_InputMap.PMap);
    
    %extract the PC list
    map_pcx(p)=MapInfo.Data_InputMap.PCX(P_r,P_c);
    map_pcy(p)=MapInfo.Data_InputMap.PCY(P_r,P_c);
    map_pcz(p)=MapInfo.Data_InputMap.DD(P_r,P_c);
    
    %extract the orientations
    map_phi1(p)=MapInfo.Data_InputMap.phi1(P_r,P_c);
    map_PHI(p)=MapInfo.Data_InputMap.PHI(P_r,P_c);
    map_phi2(p)=MapInfo.Data_InputMap.phi2(P_r,P_c);
    
    %calculate the orientation matrix
    eangs=[map_phi1(p),map_PHI(p),map_phi2(p)]*pi/180; %bruker stores in degrees
    g_bruker(:,:,p)=RTM.Rz(eangs(3))*RTM.Rx(eangs(2))*RTM.Rz(eangs(1));
    g_bruker_det(:,:,p)=g_bruker(:,:,p)*Detector_tilt; %convert to the detector
end

%% Check that this is reasonalbe

p=100;

%generate the geometry of the pattern
[EBSD_geom ] = EBSP_Gnom( RTM,[map_pcx(p) map_pcy(p) map_pcz(p)]); %you can change PC_in if you want

%update simulation for this pattern
[ Pat_sim_bruker ] = EBSP_gen( EBSD_geom,g_bruker_det(:,:,p),screen_int,screen_int.isHex ); %generate the EBSP for this iteration 

%perform a refinement
tic
[G_refined_bruker_test,~] = refine5(Pat_Ref_r(p),EBSD_geom,EBSD_geom.PC,g_bruker_det(:,:,p),SettingsXCF,screen_int,RTM_info.isHex,RTM);   
toc

%update simulation for the refined pattern
[ Pat_sim_brukerref ] = EBSP_gen( EBSD_geom,G_refined_bruker_test,screen_int,screen_int.isHex ); %generate the EBSP for this iteration 

figure;
subplot(2,2,1);
I1=pPattern(Pat_Ref_r(p).pat_in,EBSD_geom);
title('Experimental pattern');

subplot(2,2,2);
I1=pPattern(Pat_sim_bruker,EBSD_geom);
title('Simulation - Bruker measurement');

subplot(2,2,3);
I1=pPattern(Pat_sim_brukerref,EBSD_geom);
title('Simulation - Bruker measurement');

subplot(2,2,4);
I1=pPattern(Pat_sim_brukerref-Pat_sim_bruker,EBSD_geom);
title('Difference');

%% Bruker refinement of patterns

%remove the original patterns - good for memory, but plotting will need
%them to be read again..
if isfield(Pat_Ref_r,'pat_in')
    Pat_Ref_r=rmfield(Pat_Ref_r,'pat_in');
end

%preallocate variables
G_refined=zeros(3,3,numel(pat_list));
regout_Grefined=zeros(4,numel(pat_list));
disp(['Matching started']);
tic
parfor p=1:numel(pat_list)
    
    %update geometry
    [EBSD_geom ] = EBSP_Gnom( RTM,[map_pcx(p) map_pcy(p) map_pcz(p)]); %you can change PC_in if you want

    %perform a refinement
    [G_refined(:,:,p),regout_Grefined(:,p)] = refine5(Pat_Ref_r(p),EBSD_geom,EBSD_geom.PC,g_bruker_det(:,:,p),SettingsXCF,screen_int,RTM_info.isHex,RTM);
    
    %{
    if p/100 == round(p/100)
        disp(['Pattern ' int2str(p) ' run, of ' int2str(numel(pat_list))]);
    end
    %}
end
toc
disp(['All Patterns Matched']);

%% create the maps from this data

%rotate the detector data into the sample frame
G_refined_sample=zeros(3,3,numel(pat_list));

for p=1:numel(pat_list)
    G_refined_sample(:,:,p)=inv(G_refined(:,:,p)*inv(Detector_tilt));
end

%create the MTEX variables
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

prop.x = double(MapInfo.MapData.XSample); %have to transpose - not sure why...
prop.y = double(MapInfo.MapData.YSample);
prop.quality=regout_Grefined(4,:); %PH from refinement
prop.quality=prop.quality(:);

ori_Brefined = ...
  rotation('Matrix',G_refined_sample(:,:,:));
ebsd_Bruker = EBSD(ori_Brefined, ones(size(ori_Brefined)),{'notIndexed',cs_phase},prop);

figure;
colorKey = ipfHSVKey(cs_phase);

colorKey.inversePoleFigureDirection = xvector;
plot(ebsd_Bruker,colorKey.orientation2color(ebsd_Bruker.orientations))
nextAxis

colorKey.inversePoleFigureDirection = yvector;
plot(ebsd_Bruker,colorKey.orientation2color(ebsd_Bruker.orientations))
nextAxis

colorKey.inversePoleFigureDirection = zvector;
plot(ebsd_Bruker,colorKey.orientation2color(ebsd_Bruker.orientations))
nextAxis

%{
%%
%get the mean PC for the map
Mean_PC=[mean(MapInfo.Data_InputMap.PCX(:)),mean(MapInfo.Data_InputMap.PCY(:)),mean(MapInfo.Data_InputMap.DD(:))];
%generate the library for this PC
[ template_library ] = Library_Gen(EBSD_geom,screen_int,RTM_info.isHex,library_G,2,SettingsXCF);





%% Map based matching
G_SO3_map=zeros(3,3,numel(pat_list));
G_Refined_SO3_map=zeros(3,3,numel(pat_list));
G_Refined_Bruker_map=zeros(3,3,numel(pat_list));
Library_PH=zeros(numel(pat_list),1);
regout_R=zeros(4,numel(pat_list));
regout_Bruker=zeros(4,numel(pat_list));




%% Refine the Bruker map

for p=1:numel(pat_list)
    
    G_Bruker_map=gmatrix;
    PC_in=[map_pcx(p),map_pcy(p),map_pcz(p)];
    [G_Refined_Bruker_map(:,:,p),regout_Bruker(:,p)] = refine5(Pat_Ref_r(p),EBSD_geom,PC_in,gmatrix*Detector_tilt,SettingsXCF,screen_int,RTM_info.isHex,RTM);   
end



%% SO3 Search
parfor p=1:numel(pat_list)
    %compare with the library
    [G_SO3_map(:,:,p),Library_PH(p)]=fLibraryTest(template_library,library_G,Pat_Ref_r(p).FFT,SettingsXCF,XCF_data_fill,1);
   
    PC_in=[map_pcx(p),map_pcy(p),map_pcz(p)];

    %refine
    [G_Refined_SO3_map(:,:,p),regout_R(:,p)] = refine5(Pat_Ref_r(p),EBSD_geom,PC_in,G_SO3_map(:,:,p),SettingsXCF,screen_int,RTM_info.isHex,RTM);    
end

%% Image one pattern and check solutions

p=4000; %pattern number

%extract the PC for this pattern
PC_in=[map_pcx(p),map_pcy(p),map_pcz(p)];
[EBSD_geom ] = EBSP_Gnom( RTM,PC_in); %you can change PC_in if you want

%extract the pattern
Pat_EXPp=Pat_Ref_r(p).pat_in;

%create SO3 simulation
[ Pat_simp_SO3 ] = EBSP_gen( EBSD_geom,G_SO3_map(:,:,p),screen_int,screen_int.isHex ); %generate the EBSP for this iteration 
%create refined simulation
[ Pat_simp_SO3r ] = EBSP_gen( EBSD_geom,G_Refined_SO3_map(:,:,p),screen_int,screen_int.isHex ); %generate the EBSP for this iteration 

%plot the figure
f1=figure;
subplot(2,2,3);
I1=pPattern(Pat_EXPp,EBSD_geom);
title('Experimental pattern');

subplot(2,2,2);
I1=pPattern(Pat_simp_SO3-Pat_simp_SO3r,EBSD_geom);
title('Difference');

subplot(2,2,1);
I1=pPattern(Pat_simp_SO3,EBSD_geom);
title('SO(3) Template Found');

subplot(2,2,4);
I1=pPattern(Pat_simp_SO3r,EBSD_geom);
title('SO(3) Refined');
f1.Position=[100 100 800 400];

%% Convert to MTEX orientation data
G_Samp_ref=zeros(3,3,numel(pat_list));
G_Samp_SO3=zeros(3,3,numel(pat_list));

% prop.quality_so3=Library_PH; %PH from SO3
% 
%convert to the sample frame
for p=1:numel(pat_list)
    G_Samp_ref(:,:,p)=inv(G_Refined_SO3_map(:,:,p)*inv(Detector_tilt));
    G_Samp_SO3(:,:,p)=inv(G_SO3_map(:,:,p)*inv(Detector_tilt));
end

setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

% build the coordinate maps
prop.x = double(MapInfo.MapData.XSample); %have to transpose - not sure why...
prop.y = double(MapInfo.MapData.YSample);
prop.quality=regout_R(4,:); %PH from refinement
prop.quality_so3=Library_PH; %PH from SO3

ori_ref = ...
  rotation('Matrix',G_Samp_ref(:,:,:));
ori_so3 = ...
  rotation('Matrix',G_Samp_SO3(:,:,:));

ebsd_ref = EBSD(ori_ref, ones(size(ori_ref)),{'notIndexed',cs_phase},'Options',prop);
ebsd_so3 = EBSD(ori_so3, ones(size(ori_ref)),{'notIndexed',cs_phase},'Options',prop);

%% Plot some data
figure;
colorKey = ipfHSVKey(cs_phase);
colorKey.inversePoleFigureDirection = xvector;

plot(ebsd_so3,colorKey.orientation2color(ebsd_so3.orientations))
nextAxis
plot(ebsd_ref,colorKey.orientation2color(ebsd_ref.orientations))

figure;
colorKey = ipfHSVKey(cs_phase);
colorKey.inversePoleFigureDirection = yvector;

plot(ebsd_so3,colorKey.orientation2color(ebsd_so3.orientations))
nextAxis
plot(ebsd_ref,colorKey.orientation2color(ebsd_ref.orientations))

figure;
colorKey = ipfHSVKey(cs_phase);
colorKey.inversePoleFigureDirection = zvector;

plot(ebsd_so3,colorKey.orientation2color(ebsd_so3.orientations))
nextAxis
plot(ebsd_ref,colorKey.orientation2color(ebsd_ref.orientations))

figure;
plot(ebsd_so3,ebsd_so3.prop.quality_so3);
nextAxis;
plot(ebsd_ref,ebsd_so3.prop.quality);

%}