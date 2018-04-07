%script start
clear
home
close all

%load key folders into the path
Astro_FP= 'C:\Users\bbrit\Documents\MATLAB\AstroEBSD'; %please change this to your astro file path

InputUser.Mode = 'Isolated'; %modes = 
                             %Isolated = single file 
                             %Map_Single = single from map
                             %Map_All = full map, with area PC search
                             %Folder = blind folder
                             
InputUser.PatternLoc = fullfile(Astro_FP,'decks','image_Tidemo.tif');  %Input a raw image file if needed, or folder
InputUser.PatternPhase = 1; %phase number for this pattern
InputUser.PCSearch = 0;     %find the pattern centre - for single patterns
InputUser.PatternFlip = 0; %flip the pattern loaded (single & folder) - 1 = UD, 2 = LR, 3 = UD + LR

%input filename - for map related data (ignored if isolated is selected)
InputUser.EBSD_File='';
%input folder - note that pwd gives the current directory
InputUser.HDF5_folder='';
InputUser.BCF_folder='';
InputUser.EBSD_folder='';

%build the phases
InputUser.Phase_Folder = fullfile(Astro_FP,'phases');
InputUser.Phase_Input  = {'Ti_Alpha'}; %Si, Ferrite
% InputUser.Phase_Input={'Ti','Ti_Beta'};
%chose the folder for output
InputUser.FolderOut=fullfile(Astro_FP,'outputs');
InputUser.FileOut='AstroEBSD';


%PC search ranges
Settings_PCin.start=[0.5 0.35 0.6];

%% EBSD pattern settings

%background correction
Settings_Cor.gfilt=0; %use a low pass filter (do you mean high pass?)
Settings_Cor.gfilt_s=4; %low pass filter sigma

%radius mask
Settings_Cor.radius=1; %use a radius mask
Settings_Cor.radius_frac=0.95; %fraction of the pattern width to use as the mask

%hold pixel
Settings_Cor.hotpixel=0; %hot pixel correction
Settings_Cor.hot_thresh=1000; %hot pixel threshold

%reside
Settings_Cor.resize=1; %resize correction
Settings_Cor.size=150; %image width

Settings_Cor.RealBG=0; %use a real BG
Settings_Cor.EBSP_bgnum=30; %number of real pattern to use for BG

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


InputUser.DeckPath=mfilename('fullpath'); %obtain the deck path

%run the analysis code
Astro_Run;

%plot data
Astro_Plot;

%% Index the pattern
box_size=21;

peak_centre=EBSP_One.Peak_Centre;
pattern_info=EBSP_One.PatternInfo.size;

