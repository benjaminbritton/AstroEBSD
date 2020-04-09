%% Histograms of phases
% Plot histograms of cross-correlation peak heights (if required). Phase
% specific.
%TPM 07/06/19

clear xcfhist_locs xcfhist_counts xcfhist_counts_indexed xcf_locs_indexed locs quals
for iii=1:length(InputUser.Phases)
    [xcfhist_counts(iii,:),xcfhist_locs(iii,:)]=hist(RTM.Output.PeakHeight(iii,:),15);   
end

[quals,phasematch]=max(RTM.Output.PeakHeight);
for iii=1:length(InputUser.Phases)
    locs=find(phasematch==iii);
    
    if isempty(locs)
    xcfhist_counts_indexed(iii,:)=zeros(1,15);
    xcfhist_locs_indexed(iii,:)=zeros(1,15);
    continue
    end   
        
    counts=quals(locs);
    [xcfhist_counts_indexed(iii,:),xcfhist_locs_indexed(iii,:)]=hist(counts,15);
end

figure('Position', get(groot,'ScreenSize'))
subplot(1,2,1)
for iii=1:length(InputUser.Phases)
plot(xcfhist_locs(iii,:),xcfhist_counts(iii,:),'LineWidth',1.5)
hold on
end
title('All phases, all XCF values')
legend(InputUser.Phases,'Location','NorthWest')

subplot(1,2,2)
for iii=1:length(InputUser.Phases)
plot(xcfhist_locs_indexed(iii,:),xcfhist_counts_indexed(iii,:),'LineWidth',1.5)
hold on
end
title('XCF values for indexed phases')
legend(InputUser.Phases,'Location','NorthWest')

clear xcfhist_locs xcfhist_counts xcfhist_locs_indexed xcfhist_counts_indexed