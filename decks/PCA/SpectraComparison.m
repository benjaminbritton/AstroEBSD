%% Matched RC-spectra comparison to averages - for AstroEBSD v2 - T P McAuliffe 19/02/20

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clear all
home

%Choose the input patterns (RC-EBSPs) you want to compare to library.
%These will be created from the PCA_deck file
load('asdf\Results.mat');
cd('asdf\weighting1_vt0,05')

% Set up and locate plugins
InputUser.MTEX_loc='C:\Communal_MatlabPlugins\mtex-5.2.beta2';
InputUser.Astro_loc='C:\Users\tpm416\Documents\GitHub\AstroEBSD';

run([InputUser.MTEX_loc,'\startup_mtex']);
run([InputUser.Astro_loc,'\start_AstroEBSD']);

%colourmap
cmap=cbrewer('qual','Dark2',15); %can be useful to change the number in here dep on how many phases you are using

%% Average-spectra
phaselocs=cell(size(InputUser.Phases));
[~,num]=max(RTM.Output.PeakHeight,[],2);

for phase=1:size(InputUser.Phases,2)
    phaselocs{phase}=find(num==phase);
    
    meanspectra_by_phase{phase}=mean(tile.av_spectra(:,phaselocs{phase}),2);
    stdev_by_phase{phase}=std(tile.av_spectra(:,phaselocs{phase}),[],2);
    stderror_by_phase{phase}=stdev_by_phase{phase}./sqrt(length(phaselocs{phase}));

end

%
figure
for phase=1:size(InputUser.Phases,2)
    

    region(:,1)=meanspectra_by_phase{phase}-stderror_by_phase{phase};
    region(:,2)=2.*stderror_by_phase{phase};

    h=area(region,'LineStyle','none','FaceColor',cmap(phase,:));
    hold on
    h(1).FaceColor=[1,1,1];
    h(1).FaceAlpha=0;
    h(1).EdgeColor=[1,1,1];
    h(2).FaceAlpha=0.2;
    xlim([1,1000]);
   

    p(phase)=plot(meanspectra_by_phase{phase},'Color',0.8.*cmap(phase,:),'LineWidth',1);
    
    plotted(phase)=~isnan(region(1,1));

end
xlim([1,800])
ylim([0,12])

legend(p(plotted),InputUser.Phases(plotted))
print('Av_spectra','-dpng','-r300')

%% RC-Spectra
% averaging
phaselocs=cell(size(InputUser.Phases));
[~,num]=max(RTM.Output.PeakHeight,[],2);

for phase=1:size(InputUser.Phases,2)
    phaselocs{phase}=find(num==phase);
    
    meanspectra_by_phase{phase}=mean(tile.spec(:,phaselocs{phase}),2);
    stdev_by_phase{phase}=std(tile.spec(:,phaselocs{phase}),[],2);
    stderror_by_phase{phase}=stdev_by_phase{phase}./sqrt(length(phaselocs{phase}));

end

%
figure
for phase=1:size(InputUser.Phases,2)
    

    region(:,1)=meanspectra_by_phase{phase}-stderror_by_phase{phase};
    region(:,2)=2.*stderror_by_phase{phase};

    h=area(region,'LineStyle','none','FaceColor',cmap(phase,:));
    hold on
    h(1).FaceColor=[1,1,1];
    h(1).FaceAlpha=0;
    h(1).EdgeColor=[1,1,1];
    
    h(2).FaceAlpha=0.5;
    
    xlim([1,1000]);

    p(phase)=plot(meanspectra_by_phase{phase},'Color',0.8.*cmap(phase,:),'LineWidth',1);
    
    plotted(phase)=~isnan(region(1,1));

end
xlim([1,800])

legend(p(plotted),InputUser.Phases(plotted));
print('RC-spectra','-dpng','-r300')

%% RC-average comparison
% averaging
phaselocs=cell(size(InputUser.Phases));
[~,num]=max(RTM.Output.PeakHeight,[],2);

for phase=1:size(InputUser.Phases,2)
    figure
    for type=1:2
        
        if type==1
            spec=tile.spec;
            c=[1,0,0];
        else
            spec=tile.av_spectra;
            c=[0,0,1];
        end
        spec=spec./std(spec,[],1);
        
        phaselocs{phase}=find(num==phase);

        meanspectra_by_phase{phase}=mean(spec(:,phaselocs{phase}),2);
        stdev_by_phase{phase}=std(spec(:,phaselocs{phase}),[],2);
        stderror_by_phase{phase}=stdev_by_phase{phase}./sqrt(length(phaselocs{phase}));

        region(:,1)=meanspectra_by_phase{phase}-stderror_by_phase{phase};
        region(:,2)=2.*stderror_by_phase{phase};

        h=area(region,'LineStyle','none','FaceColor',c);
        hold on
        h(1).FaceColor=[1,1,1];
        h(1).FaceAlpha=0;
        h(1).EdgeColor=[1,1,1];
        
        h(2).FaceAlpha=0.2;
        xlim([1,1000]);

        p2(type)=plot(meanspectra_by_phase{phase},'LineWidth',1,'Color',0.8.*c);

        plotted(phase)=~isnan(region(1,1));
    end
    hold off

xlim([1,800])
ylim([0,12])
print([InputUser.Phases{phase},'_RC-average_comparison'],'-dpng','-r300')

legend(p2,'RCC','Average')

end

