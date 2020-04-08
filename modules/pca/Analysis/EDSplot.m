% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [mean_el,std_el]=EDSplot(InputUser,RTM,tile,quant,elements,printing,name,minV,maxV,stdev,cmap)

%%
%phase_map=zeros([size(RTM.Output.phase),length(InputUser.Phases)]);

mean_el=[];
std_el=[];
phaseid=[];
n_labels=[];

for phase_ind=1:length(InputUser.Phases)
    
phase_loc=find(RTM.Output.phase==phase_ind);

if isempty(phase_loc)
    continue
end

phase_map_sub=zeros(size(RTM.Output.phase));
phase_map_sub(phase_loc)=1;
phase_map(:,:,phase_ind)=phase_map_sub;

labels=unique(tile.map_reshaped(phase_loc));

mean_el=[mean_el, mean(quant(:,labels),2)];
std_el=[std_el, std(quant(:,labels),[],2)];
phaseid=[phaseid,phase_ind];
n_labels=[n_labels,length(labels)];
end

%%
figure('Position', get(groot,'ScreenSize'))
    
b=bar(mean_el);

for ind=1:length(phaseid)
b(ind).FaceColor=cmap(phaseid(ind),:);
b(ind).FaceAlpha=1;
hold on

set(gca,'XTickLabel',elements,'FontSize',15); 
%set(gca,'YTickLabel','')
pbaspect([2,1,1])

ctr(ind,:)=bsxfun(@plus, b(ind).XData, [b(ind).XOffset]');
y(ind,:)=b(ind).YData;


std_n(:,ind)=std_el(:,ind)./(2*sqrt(1));%n_labels(ind)));
end
er=errorbar(ctr',y',std_n);

for ind=1:length(phaseid)
er(ind).Color = 0.8*cmap(phaseid(ind),:);                            
er(ind).LineStyle = 'none';
end

legend(InputUser.Phases{phaseid})

if printing==1
savefig(['Bar_',name])
print(gcf,['Bar_',name],'-dpng','-r500');
end

%%
if isstr(minV)==1;
minV=zeros(1,length(elements));
end

if isstr(maxV)==1;
maxV=ceil((max(mean_el,[],2)))'.*ones(1,length(elements));
end

figure('Position', get(groot,'ScreenSize'))


hold on

h=radarPlot_f(1,mean_el, minV',maxV',elements,1, '-','LineWidth', 1);
if stdev==1
h1=radarPlot_f(0,mean_el+std_el./2, minV',maxV',elements,1, ':','LineWidth', 1);
h2=radarPlot_f(0,mean_el-std_el./2, minV',maxV',elements,1, ':','LineWidth', 1);
end

for ind=1:length(phaseid)
h(ind).FaceColor=cmap(phaseid(ind),:);
h(ind).FaceAlpha=0.6;
h(ind).EdgeColor='none';
h1(ind).Color=(cmap(phaseid(ind),:)).^2;
%h1(ind).EdgeColor=0.8*cmap(phaseid(ind),:);;
%h1(ind).FaceAlpha=0.3;
h2(ind).Color=(cmap(phaseid(ind),:)).^2;
%h2(ind).EdgeColor=0.8*cmap(phaseid(ind),:);;
%h2(ind).FaceAlpha=0.3;

%h(ind).FaceColor=cmap(phaseid(ind),:);
end
hold on

legend(InputUser.Phases{phaseid}, 'Location', 'Best'); 

if printing==1
savefig(['Radar_',name])
print(gcf,['Radar_',name],'-dpng','-r500');
end
end
