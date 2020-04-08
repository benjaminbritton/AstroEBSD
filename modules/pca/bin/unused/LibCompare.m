function []=LibCompare(component,EBSD_Info,RTM,InputPats)
matchedpat=zeros(EBSD_Info.PatSizeH,EBSD_Info.PatSizeW);

figure('Position', get(groot,'ScreenSize'))
subplot(2,4,1)
imagesc(InputPats(:,:,component))
title('PCA pattern')
pbaspect([1,1,1])
for iii=1:7
subplot(2,4,iii+1)
matchedpat(:,:)=RTM.matchedpats(iii,component,:,:);
imagesc(matchedpat)
title([num2str(iii), InputUser.Phases(iii),num2str(RTM.Output.PeakHeight(iii,component))])
pbaspect([1,1,1])
end
colormap('gray')

clear matchedpat component iii
end