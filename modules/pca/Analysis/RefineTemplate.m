% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [RefinedTemplatePat]=RefineTemplate(phasenum,component,RTM,tile,MapInfo,MicroscopeData,EBSD_Info,Settings_Cor,InputUser)

%run('I:\TomMcA\GitHub\RTM_indexing\start_RTM')
%run('I:\TomMcA\GitHub\AstroEBSD\start_AstroEBSD')

screensize=RTM.screensize;
Sampling_Freq=RTM.Sampling_Freq;
iterations=RTM.iterations;
LPTsize=RTM.LPTsize;
Bin_loc=RTM.Bin_loc;
Phase_Folder=RTM.Phase_Folder;

% Set up some things the same way RTM does
Data_InputMap_Start=MapInfo.Data_InputMap;
EBSD_DataInfo=MapInfo.EBSPData;
Settings_CorX=Settings_Cor;

%% Get phase number
%phasenum=1;
%component=1;

InputUser.Phase_Input=InputUser.Phases(phasenum);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( InputUser.Phase_Input,RTM.Phase_Folder);
[screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

%generate the pattern for this phase
%xi=1;yi=1;

loc=find(tile.map_reshaped==component);
loc=loc(1);
GMat_test=conv_EA_to_G([RTM.Output.euler1(loc),RTM.Output.euler2(loc),RTM.Output.euler3(loc)]);

%Define all rotation matrices needed in the code
Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

Detector_tilt = Rx(MicroscopeData.TotalTilt);
rottoplot=GMat_test*Detector_tilt;

PatternInfo.ScreenWidth=RTM.screensize;
PatternInfo.ScreenHeight=RTM.screensize;

locs=find(tile.map_reshaped==component);
PC_X = median(MapInfo.Data_InputMap.PCX(locs));
PC_Y = median(MapInfo.Data_InputMap.PCY(locs));
PC_Z = median(MapInfo.Data_InputMap.DD(locs));
pnum=median(MapInfo.Data_InputMap.PMap(locs));
PC_pat=[PC_X,PC_Y,PC_Z];

[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);

[ RefinedTemplatePat ] = EBSP_gen( EBSP_pat,rottoplot,screen_int,RTM_info.isHex ); %generate the EBSP for this iteration
end


