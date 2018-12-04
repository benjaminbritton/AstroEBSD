%script start
clear
home
close all

%load key folders into the path

%email address for status updates

InputUser.Mode = 'Map_All'; %modes = 
                             %Isolated = single file 
                             %Map_Single = single from map
                             %Map_All = full map, with area PC search
                             %Folder = blind folder
Astro_FP='C:\Users\tbritton\Documents\GitHub\AstroEBSD';                             
InputUser.PatternPhase = 1; %phase number for this pattern
InputUser.PCSearch = 1;     %find the pattern centre - for single patterns
InputUser.PatternFlip = 1; %flip the pattern loaded (single & folder) - 1 = UD, 2 = LR, 3 = UD + LR

InputUser.OnePatternPosition=[10 9]; %X and Y positions, in beam coords - for map running

%input filename - for map related data (ignored if isolated is selected)
InputUser.EBSD_File='Demo_Ben';
%input folder - note that pwd gives the current directory
InputUser.HDF5_folder='F:\Ben';
InputUser.BCF_folder='F:\Ben';

%build the phases
InputUser.Phase_Folder = [Astro_FP '\phases'];
InputUser.Phase_Input  = {'Ferrite'}; %Si, Ferrite
% InputUser.Phase_Input={'Ti','Ti_Beta'};
%chose the folder for output
InputUser.FolderOut=[Astro_FP '\outputs'];
InputUser.FileOut='AstroEBSD';


%PC search ranges
% PC_start=[0.46 0.15 0.6]; %[PCx, PCy, PCz]; %Si
Settings_PCin.start=[0.493 0.452 0.569]; %Fe
% Settings_PCin.start=[0.47 0.16 0.63]; %Co_V208C_PM\Mapping_23-04-2017
% Settings_PCin.start=[0.484,0.468,0.677];

Settings_PCin.range=[0.15 0.15 0.15]; %+- these values
Settings_PCin.array=[10 10]; %[#X,#Y] points extracted from map - will fit a PC to these points & then fit a plane

%Plotting filters - for map data
Settings_PlotFilters.MAE_Thresh=3*pi/180; %Max ok MAE, in radians
Settings_PlotFilters.IQ_Thresh= 0; %Min ok IQ

%% EBSD pattern settings

%background correction
Settings_Cor.gfilt=1; %use a low pass filter (do you mean high pass?)
Settings_Cor.gfilt_s=4; %low pass filter sigma

%radius mask
Settings_Cor.radius=1; %use a radius mask
Settings_Cor.radius_frac=0.98; %fraction of the pattern width to use as the mask

%hold pixelData_InputMap
Settings_Cor.hotpixel=0; %hot pixel correction
Settings_Cor.hot_thresh=1000; %hot pixel threshold

%resize
Settings_Cor.resize=1; %resize correction
Settings_Cor.size=150; %image width

Settings_Cor.RealBG=0; %use a real BG
Settings_Cor.EBSP_bgnum=30; %number of real pattern to use for BG

Settings_Cor.radius=1;
Settings_Cor.size=300;
Settings_Cor.SplitBG=1;
Settings_Cor.gfilt=0;
Settings_Cor.hotpixel=1;
Settings_Cor.hot_thresh=1000;

%radon searching 

%Peak Finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak hunt
Settings_Rad.max_peaks=10; %max number of peaks to return
Settings_Rad.num_peak=15; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=4; %search size in theta steps
Settings_Rad.rho_search_per=0.1; %radon search in fractions
Settings_Rad.min_peak_width=0.002; %seperation of the peak width, in pixels

%% run the code

%load key folders into the path
Init_currentdir=cd;
Init_astroloc=strfind(Init_currentdir,'AstroEBSD');
Init_path=Init_currentdir(1:Init_astroloc+8);
%check that we are in a subfolder
addpath(Init_path);
clear Init*

InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path

%run the analysis code
Astro_Run;

%plot data
Astro_Plot;

%%

cs = loadCIF('Fe-Iron-alpha.cif');
colorKey = ipfHSVKey(cs);
colorKey.inversePoleFigureDirection = xvector;
color = colorKey.orientation2color(ebsd('indexed').orientations);

% build the coordinate maps
prop.x = double(-Data_InputMap.XBeam_Map); %have to transpose - not sure why...
prop.y = double(Data_InputMap.YBeam_Map);
ori = ...
  rotation('Euler',Data_OutputMap.phi1',Data_OutputMap.PHI',Data_OutputMap.phi2');

ebsd = EBSD(ori.', ones(size(ori)),{'notIndexed',cs},'options',prop);

plot(ebsd,colorKey.orientation2color(ebsd.orientations),'micronbar','off')

saveFigure('AstroEBSD_IPF.png')

[grains,ebsd.grainId] = calcGrains(ebsd);
colorKey = axisAngleColorKey;
colorKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;

figure(1); 
color = colorKey.orientation2color(ebsd(cs.mineral).orientations);
plot(ebsd(cs.mineral),color,'micronbar','off');

hold on
plot(grains.boundary)
saveFigure('AstroEBSD_AxAngle.png')
hold off