% PCA with RTM - for AstroEBSD v2 - T P McAuliffe & T B Britton 19/02/20

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
%     MATLAB R2018a or above:
%       statistics/machine-learning, global
%       optimisation, parallel computing toolboxes.
%    
%     MTEX version 5.2.beta2 or above
%
% PCA developed by McAuliffe et al [1] following Wilkinson et al [2].
% RTM developed by Foden et al [3].

% Example DATA is available at https://doi.org/10.5281/zenodo.3737987.

%% clear and reset
clear variables; close all;

%% DIRECTORIES info %%
InputUser.HDF5_folder='E:\Tom'; %folder where the h5/bcf is located

% Cell array of .H5 files to run PCA on: can be everything in the
% directory, or specific filenames. Example syntax:
files= {'SuperalloyExample.h5'};

% Can also do:
%filelist=dir([InputUser.HDF5_folder,'/*.h5']);
%files = {filelist.name}.';
%clear filelist

% This is the folder in which the outputs will be saved.
InputUser.SavingFolder='C:\Users\tpm416\Documents\PCA_Output';

% Plugin locations
InputUser.MTEX_loc='C:\Communal_MatlabPlugins\mtex-5.2.beta2';
InputUser.Astro_loc='C:\Users\tpm416\Documents\GitHub\AstroEBSD';

%% CRITICAL settings and setup - these determine which analysis to run.
% YES = 1, NO = 0

%These need to be the names of '.pha' files in the '\phases\phasefiles'
%folder. There need to be corresponding '.cif' and '.BIN' files in 'cifs'
%and 'masterpatterns' folders. - Only relevant to RTM analyses.
InputUser.Phases={'Ni','ZrC','M6C'};

% Run PCA on EBSD, EDX, or both
PCA_Setup.PCA_EBSD=1;
PCA_Setup.PCA_EDX=1;

% Do it with a relative weighting (applied to the EBSD data)?
PCA_Setup.weighted=1;
PCA_Setup.EBSD_weighting=1; %EBSD weighting. Can give as a vector and script will loop. (1)

% Run the RTM (requires EBSD PCA)
PCA_Setup.RTM=1;

% Introduce a spatial weighting kernel after Guo et al [4]
% This will be convolved with the data matrix with stride 1 - each pattern
% is a weighted sum of its kernel. Kernel fn can be adjusted in LL settings.
PCA_Setup.SpatialKernel=1;
PCA_Setup.KernelRadius=3;%n in pixels - must be ODD

% Pattern centre refinement settings
Refine.run=1; %Do you want RTM PC refinement to run? If = 0 uses .h5 loaded PC model.
Refine.phase=1; %Which phase number do you want to select points for? (Number in InputUser.Phases)

%% MEDIUM LEVEL settings - you may fairly frequently want to change these.
printing=1; % Do you want files to display and save?

% Dimension of square tiling (3 = 3x3 grid) %select '1' if no tiling wanted
% This has a big effect on RAM
PCA_Setup.crop_factor=3;

% This systematically chooses number of components to retain in the tiles
PCA_Setup.variance_tolerance=0.1;% in percent. (0.1)

% Square dimension EBSPs are cropped to - this has a big effect on RAM
RTM_setup.screensize=200; %(120 - 300 usual)

%% LOWEST LEVEL settings - most users will not need to adjust these.

% Extra PCA options
PCA_Setup.MedFilt=0;% Median filter label map (0)
PCA_Setup.components=[];% if not blank, then will override number of components to retain with this number (for every tile).
PCA_Setup.run_one_tile=[]; %if not blank, pipeline will ONLY run on this tile.

% RTM procedure settings
RTM_setup.Sampling_Freq=8; %of SO3 space (7 - 10)
RTM_setup.iterations=4; %(4)
RTM_setup.LPTsize=500; %(500)
RTM_setup.parsearch=1; %0: no parallelisation ; 1: search template library in parallel ; 2: search patterns in parallel

% Spatial kernel settings
PCA_Setup.KernelFnLabel='Guo'; %for labelling
PCA_Setup.KernelFunction = @(distance,r) ((1 - (distance)./r).^2).^2;
%check with k=GenerateKernel(PCA_Setup.KernelFunction,PCA_Setup.KernelRadius);

% AstroEBSD based background correction
Settings_Cor.gfilt=1; %use a high pass filter (1)
Settings_Cor.gfilt_s=4; %low pass filter sigma (4)
Settings_Cor.radius=1; %use a radius mask (0)
Settings_Cor.radius_frac=0.85; %fraction of the pattern width to use as the mask (0.85)
Settings_Cor.hotpixel=0; %hot pixel correction (0)
Settings_Cor.hot_thresh=1000; %hot pixel threshold (1000)
Settings_Cor.resize=1; %resize correction (1)
Settings_Cor.size=RTM_setup.screensize; %image height
Settings_Cor.SplitBG=1; %split chip fwix (1)
Settings_Cor.RealBG=0; %use a real BG (0)
%Settings_Cor.EBSP_bgnum=30; (30)
Settings_Cor.Square=1; %square crop (1)
Settings_Cor.SquareCrop=1; %(1)
Settings_Cor.LineError=1; %(1)
Settings_Cor.MeanCentre=1; %(1)
Settings_Cor.channum=2048; % channels in the EDS spectra (if any) - usually 2048

% PC Refinement specifics
Refine.num_pts=6; %number of points to use for the refinement step
Refine.ss=0.08; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=10; %number of interations
Refine.debug=1; %plot the EBSPs if you need to debug

%% Add packages to the path
run([InputUser.MTEX_loc,'\startup_mtex']);
run([InputUser.Astro_loc,'\start_AstroEBSD']);

%% Start the anaysis

%... sanity check your settings... this is helpful but won't catch everything!
SettingsReadout
%... and run
PCArun_withRTM_seriesloading 

%% References
% [1] doi.org/10.1016/j.ultramic.2020.112944
% [2] doi.org/10.1016/j.ultramic.2018.09.011
% [3] doi.org/10.1016/j.ultramic.2019.112845
% [4] doi.org/10.1080/10618600.2014.912135

