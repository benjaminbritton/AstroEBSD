folders=["I:\TomMcA\PCA_output\forpaper_2"];
cd(folders)
list=dir('*_EBSD*');
addpath(genpath('I:\TomMcA\GitHub'),'-end')
ind_forbidden=[12,19];

name={};
comps=[];
skip=0;
inds=[];
phases=[];

mean_PHs=[];
std_PHs=[];
mean_ers=[];

%%
for ind=1:3%length(list)
    
    if ind == 10 | ind == 12 | ind == 29
        skip=skip+1;
        continue
        
    end
        
cd(folders)
cd(list(ind).name)
cd('weighting1_vt0,1')

clearvars -except backfol homefolder list folders ind comps name skip inds phases mean_PHs std_PHs mean_ers

phasepresence=[];

%Load data
load('Results.mat')
%pRTM(RTM,tile,EDSWindowing,1,InputUser)
%print(gcf,['Output',num2str(ind)],'-dpng','-r500');
close all
homefolder=pwd;

for p=1:13
    if size(find(RTM.Output.phase==p),1)>0;
        phasepresence(p)=1;
    else 
        phasepresence(p)=0;
    end
end

% cd([homefolder,'\AverageSpectra'])
% % loads a file called 'results'
% [quantif,elements]=Spectra_PCA(pwd);

%assemble elements in correct order
% for elcheck=1:13 
%     els_tocheck = {'B','C','O','Al','S','Cr','Co','Ni','Zr','Mo','Ta','W'};
%     if isequal(elements(elcheck),els_tocheck(elcheck))==0
%         
%         if elcheck>1
%         quantif=[quantif(:,1:elcheck-1),zeros(size(quantif,1),1),quantif(:,elcheck:end)];
%         elements={elements{1:elcheck-1},els_tocheck{elcheck},elements{elcheck:end}};
% 
%         else 
%             
%         quantif=[zeros(size(quantif,1),1),quantif(:,elcheck:end)];
%         elements={els_tocheck{elcheck},elements{elcheck},elements{elcheck+1:end}};
% 
%         end
%         
%     end 
% end
% elements=els_tocheck;

% if isequal(elements(6),{'Ti'})
%     
% else
%     elements=[elements(1:5),'Ti',elements(6:end)];
%     quantif=[quantif(:,1:5),zeros(size(quantif,1),1),quantif(:,6:end)];
% end

% subtracted=0;
% for index=1:length(elements)
%     ind_cols=index-subtracted;
% if quantif(:,ind_cols)==zeros(size(quantif,1),1)
%     quantif(:,ind_cols)=[];
%     elements(ind_cols)=[];
%     subtracted=subtracted+1;
% end
% end
% clear subtracted

% quantif=quantif';
mean_el=[];
std_el=[];
phaseid=[];
n_labels=[];
mean_PH=zeros(13,13,1);
std_PH=zeros(13,13,1);
mean_er=zeros(13,13,1);


for phase_ind=1:length(InputUser.Phases)
    
phase_loc=find(RTM.Output.phase==phase_ind);

% if isempty(phase_loc)
%     mean_el=[mean_el,zeros(13,1)];
%     continue
% end

phase_map_sub=zeros(size(RTM.Output.phase));
phase_map_sub(phase_loc)=1;
phase_map(:,:,phase_ind)=phase_map_sub;

labels=unique(tile.map_reshaped(phase_loc));

% mean_el=[mean_el, mean(quantif(:,labels),2)];
% std_el=[std_el, std(quantif(:,labels),[],2)];
phaseid=[phaseid,phase_ind];
n_labels=[n_labels,length(labels)];

%Get RTM peak height info
mean_PH(:,phase_ind)=mean(RTM.Output.PeakHeight(:,labels),2);
std_PH(:,phase_ind)=std(RTM.Output.PeakHeight(:,labels),[],2);
mean_er(:,phase_ind)=std_PH(:,phase_ind)./sqrt(length(labels));

end

comps=cat(3,comps,mean_el);
name=[name,InputUser.HDF5_file];
inds=[inds,ind];
phases=cat(1,phases,phasepresence);

mean_PHs=cat(3,mean_PHs,mean_PH);
std_PHs=cat(3,std_PHs,std_PH);
mean_ers=cat(3,mean_ers,mean_er);

%clearvars -except backfol homefolder list folders ind comps name skip inds phases mean_PHs std_PHs mean_ers labels
end