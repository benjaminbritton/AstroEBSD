%% PCA postanalysis code  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

if exist('inds')
    list_to_run=length(inds);
else
    list_to_run=length(list);
end
    
for ind=1:(list_to_run)
    
disp(['Analysing ',num2str(ind), ' of ',num2str((list_to_run))]);    

if exist('folders')
    cd(folders)
    cd(list(ind).name)
    cd(namesetting)

    if ismember(ind,ind_forbidden)==1
        continue
    end
    
else
    cd(list{ind})
    cd(namesetting)
end

clearvars -except backfol singletile homefolder list folders ind ind_forbidden RTM_update post cmap namesetting inds InputUser_update list_to_run InputUser_update

%Load data
load('Results.mat')
%pRTM(RTM,tile,EDSWindowing,1,InputUser)
%print(gcf,['Output',num2str(ind)],'-dpng','-r500');
close all
homefolder=pwd;

% %update bin_locs?
% if RTM_update.updatelocs==1
%     %fix for having changed phase folders
%     RTM.Bin_loc=RTM_update.Bin_loc;
%     RTM.Phase_Folder=RTM_update.Phase_folder;
% end
% 
% %update data location?
% if InputUser_update.update==1
%     %fix for having changed data folders
%     InputUser.HDF5_folder=InputUser_update.HDF5_folder;
% end

if exist('PCA_Setup')==1 & length(PCA_Setup.components)==1
    post.singletile=1;
end

if post.singletile==1
    tile.xtile=[1];
end

%Generate analysis folder in saving folder
analysisfolder=fullfile(homefolder,'Analysis');
mkdir(analysisfolder)
cd(analysisfolder)

%% Create EBSD plots
if post.EBSD==1
    mkdir('IPFMaps')
    cd('IPFMaps')
    ntiles=length(tile.xtile);
    [OutputUser]=RTM_save(RTM,ntiles,InputUser,MapInfo,MicroscopeData);
    [EBSD_template,cs2]=RTM_open(OutputUser,RTM);
    Astro_IPF_plot(EBSD_template,1)
    
    post.EBSD=EBSD_template;
    
    cd(analysisfolder)
    close all
end

%cmap(6,:)=cmap(4,:);
%cmap(7,:)=[0.8,0,0];

%% Bruker element quantification of RC spectra
if post.RCSpectra==1
    cd([homefolder,'\RCSpectra'])
    % loads a file called 'results'
    [quant,elements]=Spectra_PCA(pwd); 
    subtracted=0;
    cd(analysisfolder)

    %remove zeros as this causes problems with the radar plot
%     for index=1:length(elements)
%         ind_cols=index-subtracted;
%     if quant(:,ind_cols)==zeros(size(quant,1),1)
%         quant(:,ind_cols)=[];
%         elements(ind_cols)=[];
%         subtracted=subtracted+1;
%     end
%     end
%     clear subtracted

    quant(quant==0)=1e-10;

    mkdir('BrukerQuant')
    cd('BrukerQuant')
    %elements in columns
    %minV='mins';
    %maxV='maxs';

    % maxima and minima for radar plot. If want auto scaling any string is
    % fine.
    %minV=[0,0,0,0,0,0,0,0,0,0,0,0];
    %maxV=[55,2,10,2,15,30,30,5,5,15,10,5];
    
    minV=post.chemlims.minV;
    maxV=post.chemlims.maxV;

    stdev=post.RCSpectra_std;

    [mean_el,std_el]=EDSplot(InputUser,RTM,tile,quant',elements,1,['BrukerQuant_RC_',num2str(ind)],minV,maxV,stdev,cmap);
    
    post.mean_el_rc=mean_el;
    post.std_el_rc=std_el;
    
    close all
    ElementMaps(tile,quant',elements,1)
    close all
    cd(analysisfolder)
end

%% Bruker element quantification of average spectra
if post.AvSpectra==1
    cd([homefolder,'\AverageSpectra'])
    % loads a file called 'results'
    [quant,elements]=Spectra_PCA(pwd); 
    subtracted=0;
    cd(analysisfolder)
    
    %remove zeros as this causes problems with the radar plot 
%     for index=1:length(elements)
%         ind_cols=index-subtracted;
%     if quant(:,ind_cols)==zeros(size(quant,1),1)
%         quant(:,ind_cols)=[];
%         elements(ind_cols)=[];
%         subtracted=subtracted+1;
%     end
%     end
%     clear subtracted

    quant(quant==0)=1e-10;

    mkdir('BrukerQuant_Average')
    cd('BrukerQuant_Average')
    %elements in columns
    %minV='mins';
    %maxV='maxs';

    % maxima and minima for radar plot. If want auto scaling any string is
    % fine.
    %minV=[0,0,0,0,0,0,0,0,0,0,0,11];
    %maxV=[55,2,10,2,15,30,30,5,5,15,10,11];
    minV=post.chemlims.minV;
    maxV=post.chemlims.maxV;
    
    comb=post.comb;

    stdev=post.AvSpectra_std; % wheteher to plot the standard deviations or not
    [mean_el,std_el]=EDSplot(InputUser,RTM,tile,quant',elements,1,['BrukerQuant_Average_',num2str(ind)],minV,maxV,stdev,cmap);
    
    post.mean_el_av=mean_el;
    post.std_el_av=std_el;
   
    
    close all
    ElementMaps(tile,quant',elements,1)
    ElementMaps_combined(tile,quant',elements,comb,1);
    close all
    cd(analysisfolder)
end

if post.phase==1
    phase=RTM.Output.phase;
    phase(phase==0)=nan;
    
    imagesc(phase,'AlphaData',~isnan(phase))
    colormap(gca,cmap)
    title('RTM phase assignment')
    hcb1=colorbar;
    pbaspect([1.33,1,1])
    set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
    set(gca,'XTick',''); set(gca,'YTick','')
    caxis(gca,[1,length(InputUser.Phases)]);
    savefig('Phases')
    print(gcf,'Phases','-dpng','-r300')
end



end
%%
clearvars -except backfol homefolder list folders ind ind_forbidden RTM_update post cmap namesetting inds



            
    

