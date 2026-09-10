%% AstroECP v1
% https://arxiv.org/abs/2507.00354

% -------
% a digital twin of our SEM within a graphical user interface to help navigate the electron channeling patterns 
% and quickly access specific ECCI conditions for subsequent analysis, based on the approaches and functions used in 
% AstroEBSD (Britton, Tong et al., 2018) and MTEX v5.11.1 (Bachmann et al., 2010; Hielscher et al., 2019). 
% AstroECP reads the dynamically calculated Kikuchi reference patterns generated from pattern simulation software 
% that include: MapSweeper EBSD dynamical patterns (higher quality Bloch Wave Kikuchi Diffraction (BWKD) patterns, 
% the open-source software EMsoft (Singh, Marc de Graef et al., 2017) and also Bruker DynamicS. The GUI then 
% reprojects the channeling pattern based upon the crystal orientation and ECP conditions. Manipulation of the 
% simulation can be performed through virtual movement of the stage and microscope to inform movement of the sample 
% within the SEM. This helps to effectively orient the crystal for the development of crystallographic contrast with 
% a specific [uvw] direction along the optic axis of the SEM, and provide imaging conditions that are suitable for 
% the collection of ECCI micrographs.
% -------- 

clear; home; close all;

%% set up the file locations 
Input_Data.astro_location='C:\Users\benja\OneDrive\Documents\GitHub\AstroEBSD\';
Input_Data.mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.11.1\'; 

%% ECP location and information

%%  Tescan Data - comment out this block if you want to run TFS data

Input_Data.image_folder=[Input_Data.astro_location '\modules\AstroECP'];
Input_Data.image_name='Si_SAECP_example.tif';     %available in https://github.com/ExpMicroMech/AstroEBSD/blob/main/modules/AstroECP/
Input_Data.image_frame=1;     %frame number for Tescan data, this variable should not exist for other file types
Input_Data.ECP_type='TESCAN'; %supported types: 'TESCAN', 'TFS', 'other'

%values for the example data
Input_Data.PC_in=[0.5 0.5 3.9031]; % starting PC - AstroEBSD convention [PCx, PCy, DD]
Input_Data.eangs=[87.4431,0.674269,-96.3694]; % for the example pattern 'Si_SAECP_example.tif'


%set the greyscale colorlim for the ECP
% Input_Data.ECP_Pat_clim=[2 5]; % default settings of histogram - regular patterns
Input_Data.ECP_Pat_clim=[-0.5 1]; % default settings of histogram - used for the BG corrected one

%% TFS data - comment this block out if you want to run the Tescan data

% Input_Data.image_folder=[Input_Data.astro_location '\modules\AstroECP'];
% Input_Data.image_name='20keV_3.2nA_Ref.tif';     %available in https://github.com/ExpMicroMech/AstroEBSD/blob/main/modules/AstroECP/
% Input_Data.ECP_type='TFS'; %supported types: 'TESCAN', 'TFS', 'other'
% 
% %values for the example data
% Input_Data.PC_in=[0.5 0.5 8.39949];
% Input_Data.eangs=[-92.4754,1.0886,133.8429]; % for the example pattern '20keV_3.2nA_Ref.tif'
% 
% %set the greyscale colorlim for the ECP
% Input_Data.ECP_Pat_clim=[1 4]; % default settings of histogram

%% Other ECP information 

Input_Data.Stage_in=[0 0 0]; % stage rotation settings, in degrees [Rx, Ry, Rz]


%% Crystal information

Input_Data.Phase_Input{1}='Si_20kV'; %pick the phase - this should be in \AstroEBSD\phases\phasefiles

Input_Data.crystal_shape='cube'; %the crystal shape to use from MTEX - 'hex' and 'cube' coded or 'orthorhombic'
Input_Data.miller1=[0 0 1]; %miller indicies for PF plot of the crystal
Input_Data.miller2=[1 1 0]; %miller indicies for PF plot of the crystal

%location of the phase information
Input_Data.Phase_Folder = fullfile(Input_Data.astro_location,'phases'); %location of the AstroEBSD phases super-folder
Input_Data.Bin_loc = fullfile(Input_Data.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM

Input_Data.Index=1; %show the Atlas and Index buttons
%% Start MTEX and AstroEBSD
try EBSD;
catch
    run(fullfile(Input_Data.mtex_location,"startup_mtex.m"));
end

try astro_loadcheck
catch
    run(fullfile(Input_Data.astro_location,"start_AstroEBSD.m"));
end

%% set some MTEX preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');

%% Load data and run the GUI

% Experimental ECP pattern load
[ECP_Pat] = AstroECP_Load(Input_Data);

% Now Run the AstroECP GUI
[Output_Data]=fAstroECP(Input_Data,ECP_Pat);

