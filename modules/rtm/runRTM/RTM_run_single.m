% This code is copyright Alex Foden and Ben Britton 09/04/2019
% Do not distribute.
%
% This RTM-EBSD (refined template indexing) code and its associated scripts
% may only be shared with express and direct permission of
% Alex Foden & Ben Britton
% 
% If you would like a copy or access to the repository
% Please contact b.britton@imperial.ac.uk or a.foden16@imperial.ac.uk
%
% Ben Britton has a list of authorised users of this code
% 
% Do not fork
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

%% Establish key variables & anon functions for the code to run
% InputUser.BCF_folder=InputUser.HDF5_folder; %this is used for the BCF reader

%build the phases using Refine Template Matching - this is a modified version of the AstroEBSD
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,Phase_Folder, Bin_loc );

%load the crystal structure from the CIF file - this will be used for the symmetry elements
cs=loadCIF(RTM_info.cif_file); %CIF file for the cyrstal you aare indexing. This is for symmetry in generation of the library

%Define all rotation matrices needed in the code
Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%Set up the FFT filters for FFT and LPT
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( screensize, LPTsize );

%Create the master pattern used to generate the library
[screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);


%% load the H5 file
%Unpack the data
% pTime('Reading HDF5 / BCF data',time1);
% [MapData,MicroscopeData,Phase_BrukerData,EBSD_DataInfo] = BCF_HDF5( InputUser );
% [Data_InputMap_Start] = EBSD_Map(MapData,MicroscopeData);

% pTime('Updating the PC map to adjust for square data',time1);
% %update the pattern centre for square cropping
% [ Data_InputMap ] = PC_square( EBSD_DataInfo, Data_InputMap_Start,Settings_CorX );

%Define the tilt of the detector
Detector_tilt = Rx(Detecrtor_tilt);

%% build the SO(3) library

PatternInfo.ScreenWidth=screensize;
PatternInfo.ScreenHeight=screensize;

%Set up the screen
[ EBSP_av ] = EBSP_Gnom( PatternInfo,PC );

pTime('Template library generation',time1);

%Generate the library orientations
[ library_G ] = SO3_rotmat_gen( cs,Sampling_Freq);
pTime('orintations generated',time1);
% generate the templates
[ template_library ] = Library_Gen(EBSP_av,screen_int,RTM_info.isHex,library_G,XCF_type,SettingsXCF); %Need to update to get rid of XCF_type

pTime('library complete',time1);
%% Now match the patterns


pTime('Matching pattern',time1);


    
    %build the geometry
    [ PatternCor,PatternInfo ] = EBSP_BGCor( Input_EBSP,Settings_CorX );
    [ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC);
    
    %match from library
    [rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, PatternCor, library_G, template_library); %update to remove shifts, we don't need them now
    
    %refine the orientation
    [rotmat_best,regout] = refine4(PatternCor,PC,EBSP_pat,rotmat_1st,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int,RTM_info.isHex,LPTsize,iterations); 
    
    %convert to euler angles
    eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
    TemData=[eang_template regout]';
    
