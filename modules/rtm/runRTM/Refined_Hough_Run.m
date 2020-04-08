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
InputUser.BCF_folder=InputUser.HDF5_folder; %this is used for the BCF reader

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
pTime('Reading HDF5 / BCF data',time1);
[MapData,MicroscopeData,Phase_BrukerData,EBSD_DataInfo] = BCF_HDF5( InputUser );
[Data_InputMap_Start] = EBSD_Map(MapData,MicroscopeData);

pTime('Updating the PC map to adjust for square data',time1);
%update the pattern centre for square cropping
[ Data_InputMap ] = PC_square( EBSD_DataInfo, Data_InputMap_Start,Settings_CorX );

%Define the tilt of the detector
Detector_tilt = Rx(MicroscopeData.TotalTilt);

%% Now refine the orientations
TemData=zeros(7,Data_InputMap.max_pats);

pTime('Matching patterns',time1);
parfor pat=1:Data_InputMap.max_pats
% for pat=1:10
    %Extract the pattern to index
    pnum=Data_InputMap.PMap(pat);
    EBSP_raw = bReadEBSP(EBSD_DataInfo,pnum);
    
    %use the pattern centre model
    PC_X = Data_InputMap.PCX(pat);
    PC_Y = Data_InputMap.PCY(pat);
    PC_Z = Data_InputMap.DD(pat);
    
    PC_pat=[PC_X,PC_Y,PC_Z];
    
    %build the geometry
    [ PatternCor,PatternInfo ] = EBSP_BGCor( EBSP_raw,Settings_CorX );
    [ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);
    
    %Get the Hough rottation information
    phi1 = Data_InputMap.phi1(pat);
    PHI = Data_InputMap.PHI(pat);
    phi2 = Data_InputMap.phi2(pat);
    
    rotmat_1st = Rz(phi2 * pi/180)*Rx(PHI * pi/180)*Rz(phi1 * pi/180)*Detector_tilt;
    %refine the orientation
    [rotmat_best,regout] = refine4(PatternCor,PC_pat,EBSP_pat,rotmat_1st,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int,RTM_info.isHex,LPTsize,iterations); 
    
    %convert to euler angles
    eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
    TemData(:,pat)=[eang_template regout]';
    
    if pat/100 == 0
    disp(['pattern ' num2str(pat) ' of ' num2str(Data_InputMap.max_pats)])
    end
    
end