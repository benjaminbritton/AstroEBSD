% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
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

% Updated 13/02/20 TPM to include in AstroEBSD version 2.0

clear % clear any previous data
home % clear the command window
close all %Close all open figures
time1=clock; %create a clock for pTime to work

%% setup

InputUser.MTEX_loc='C:\Communal_MatlabPlugins\mtex-5.2.beta2';
InputUser.Astro_loc='C:\Users\tpm416\Documents\GitHub\AstroEBSD';

%% Add MTEX and AstroEBSD to the path

run(fullfile(InputUser.MTEX_loc,'startup_mtex.m')); %start MTEX
run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m')); %start astro

%% Inputs

InputUser.HDF5_folder='E:\Tom'; %folder location for the BCF/HDF5 file
InputUser.HDF5_file = 'SuperalloyExample.h5'; %name of the BCF/HDF5 file
InputUser.Phases  = {'Ni','ZrC'}; %name of the phase you are indexing

RTM.screensize = 128; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=7; %Set the SO(3) sampling freq in degrees
RTM.iterations = 4;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 500; %LPT size used in pixels

%% Set up background corrections for reference image

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

%% now run the code
RTM.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'masterpatterns'); %location of the binary files used for RTM
RTM.parsearch=2;

[MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadHDF5( InputUser );
[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);

%% Run the matching

pTime('Starting RTM matching',time1);
RTM=RTM_run(InputUser,RTM,Settings_CorX,MicroscopeData,MapInfo,time1);

%% Assign the phases and quality
[RTM.Output.qual,RTM.Output.phase]=max(RTM.Output.PeakHeight,[],2);

for i=1:MicroscopeData.NPoints
    RTM.Output.euler1(i,:)=RTM.Output.Eulers(1,i,RTM.Output.phase(i));
    RTM.Output.euler2(i,:)=RTM.Output.Eulers(2,i,RTM.Output.phase(i));
    RTM.Output.euler3(i,:)=RTM.Output.Eulers(3,i,RTM.Output.phase(i));
end

%% Save the data
pTime('Saving RTM data',time1);

InputUser.EBSD_File=InputUser.HDF5_file;
MapInfo.cropped_height=MicroscopeData.NROWS;
MapInfo.cropped_width=MicroscopeData.NCOLS;

OutputUser=RTM_save(RTM,1,InputUser,MapInfo,MicroscopeData);

OutputUser.Phases=InputUser.Phases;
[EBSD_template,~]=RTM_open(OutputUser,RTM);

printing=0;
Astro_IPF_plot(EBSD_template,printing);
