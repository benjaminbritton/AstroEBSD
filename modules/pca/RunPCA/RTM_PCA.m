%% Run the RTM on RC-EBSPs  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clear testArray

%Take input patterns from tile structure (RC-EBSPs)
RTM.InputPats=tile.pat;
ComponentMap=tile.map_reshaped;

%This is used for

% if PCA_Setup.PCA_EDX==1
%     ElementMaps=EDSWindowing;
% end
%     
%% Run the RTM
%matchedpats=zeros(length(InputUser.Phases),size(InputPats,3),EBSD_Info.PatSizeH,EBSD_Info.PatSizeW);

% Template match
RTM.ComponentMap=tile.map_reshaped;

RTM=RTM_run(InputUser,RTM,Settings_Cor,MicroscopeData,MapInfo,t1); %input patterns exist so should default to using these.
    
%% Assign to PCA label locations etc
[phaseassign_val,phaseassign_loc]=max(RTM.Output.PeakHeight,[],2);
total_ret_comps=sum(tile.ret_comps);

for a=1:total_ret_comps
bestphase=phaseassign_loc(a);
Euler1_assign(a)=RTM.Output.Eulers(1,a,bestphase);
Euler2_assign(a)=RTM.Output.Eulers(2,a,bestphase);
Euler3_assign(a)=RTM.Output.Eulers(3,a,bestphase);
end

RTM.Output.qual=zeros(size(ComponentMap));
for a=1:total_ret_comps
    locs=find(ComponentMap==a);
    RTM.Output.qual(locs)=phaseassign_val(a);
    clear locs
end

RTM.Output.phase=zeros(size(ComponentMap));
for a=1:total_ret_comps
    locs=find(ComponentMap==a);
    RTM.Output.phase(locs)=phaseassign_loc(a);
    clear locs
end

RTM.Output.euler1=zeros(size(ComponentMap));
for a=1:total_ret_comps
    locs=find(ComponentMap==a);
    RTM.Output.euler1(locs)=Euler1_assign(a);
    clear locs
end

RTM.Output.euler2=zeros(size(ComponentMap));
for a=1:total_ret_comps
    locs=find(ComponentMap==a);
    RTM.Output.euler2(locs)=Euler2_assign(a);
    clear locs
end

RTM.Output.euler3=zeros(size(ComponentMap));
for a=1:total_ret_comps
    locs=find(ComponentMap==a);
    RTM.Output.euler3(locs)=Euler3_assign(a);
    clear locs
end

clear loc Euler1_assign Euler2_assign Euler3_assign phaseassign_loc phaseassign_val a total_ret_comps ComponentMap

%RTM.Output=RTMOutput;
%RTM.matchedpats=matchedpats;

%% Print results from the RTM
if PCA_Setup.PCA_EDX==1
    pRTM(RTM.Output,tile,[],printing,InputUser,PCA_Setup)
else
    pRTM(RTM.Output,tile,[],printing,InputUser,PCA_Setup)
end
close all
