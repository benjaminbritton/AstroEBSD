% created by Alex Foden - modified T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [EBSD_template,cs2]=RTM_open(OutputUser,RTM)

OutputUser.Phase_Input=OutputUser.Phases;

[ RTM_MapData,RTM_MicroscopeData,RTM_PhaseData,RTM_EBSPData ] = bReadHDF5( OutputUser );

for i =1:length(OutputUser.Phases)
    [~,~,~,~,~,RTM_newinfo_it ] = Phase_Builder_RTM( OutputUser.Phases(i),RTM.Phase_Folder );
    RTM_newinfo.cif_file{i}=RTM_newinfo_it.cif_file;
    RTM_newinfo.bin_file{i}=RTM_newinfo_it.bin_file;
    RTM_newinfo.isHex{i}=RTM_newinfo_it.isHex;
end

r_pl=reshape(RTM.Output.phase,[],1); %reshaped phase list

propsT.x=RTM_MapData.XSample;
propsT.y=RTM_MapData.YSample;
propsT.quality=RTM_MapData.PeakHeight;
phasesT=propsT.quality*0;

rot_Template=rotation('Euler',RTM_MapData.phi1*degree,RTM_MapData.PHI*degree,RTM_MapData.phi2*degree,'ZXZ');

%%
for phaselist=1:length(OutputUser.Phases)
cs=loadCIF(RTM_newinfo.cif_file{phaselist}); %CIF file for the cyrstal you aare indexing. This is for symmetry in generation of the library
cs.mineral=OutputUser.Phases{phaselist};

cs2{phaselist}=cs;
end

EBSD_template=EBSD(rot_Template,r_pl,cs2,'options',propsT);


end