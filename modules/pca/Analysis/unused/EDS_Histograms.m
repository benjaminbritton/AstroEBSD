function [] = EDS_Histograms(InputUser,RTM,tile,EDSWindowing,printing,fileno)

%%
phase_map=zeros([size(RTM.Output.phase),length(InputUser.Phases)]);
cmap=cbrewer('qual','Dark2',10);
mean_pkh=[];
std_pkh=[];
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

mean_pkh=[mean_pkh, mean(EDSWindowing.p(:,labels),2)];
std_pkh=[std_pkh, std(EDSWindowing.p(:,labels),[],2)];
phaseid=[phaseid,phase_ind];
n_labels=[n_labels,length(labels)];
end

%%
figure('Position', get(groot,'ScreenSize'))
    
b=bar(mean_pkh);

for ind=1:length(phaseid)
b(ind).FaceColor=cmap(phaseid(ind),:);
b(ind).FaceAlpha=1;
hold on

set(gca,'XTickLabel',EDSWindowing.elements,'FontSize',15); 
%set(gca,'YTickLabel','')
pbaspect([2,1,1])

ctr(ind,:)=bsxfun(@plus, b(ind).XData, [b(ind).XOffset]');
y(ind,:)=b(ind).YData;


std_pkh_n(:,ind)=std_pkh(:,ind)./(2*sqrt(1));%n_labels(ind)));
end
er=errorbar(ctr',y',std_pkh_n);

for ind=1:length(phaseid)
er(ind).Color = 0.8*cmap(phaseid(ind),:);                            
er(ind).LineStyle = 'none';
end

legend(InputUser.Phases{phaseid})

if printing==1
print(gcf,['EDS_Hist',num2str(fileno)],'-dpng','-r500');
end

%%
minV=zeros(1,14);
maxV=ceil((max(mean_pkh,[],2)))'.*ones(1,14);
figure('Position', get(groot,'ScreenSize'))


hold on

h=radarPlot_f(1,mean_pkh, minV',maxV',EDSWindowing.elements,1, '-','LineWidth', 1);
h1=radarPlot_f(0,mean_pkh+std_pkh./2, minV',maxV',EDSWindowing.elements,1, ':','LineWidth', 1);
h2=radarPlot_f(0,mean_pkh-std_pkh./2, minV',maxV',EDSWindowing.elements,1, ':','LineWidth', 1);

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
print(gcf,['EDS_Radar',num2str(fileno)],'-dpng','-r500');
end


end