%This plots the full AOI and compares to reassembled watershed selections
%TPM 07/06/19

figure('Position', get(groot,'ScreenSize'))
cmap=cbrewer('qual','Set1',100);
%figure('Position',[100,350,1500,600])
[hsub,wsub]=tight_subplot(1,4,0.03,0.03);

axes(hsub(1))
imagesc(Plots.Radon)
colormap(gca,'gray')
%hcb1=colorbar;
pbaspect([1.333,1,1])

axes(hsub(2))
imagesc(tile.watershed_reshaped)
colormap(gca,cmap)
%hcb2=colorbar;
pbaspect([1.333,1,1])

axes(hsub(3))
imagesc(tile.map_reshaped_original)
colormap(gca,cmap)
%hcb2=colorbar;
pbaspect([1.333,1,1])

axes(hsub(4));
imagesc(tile.map_reshaped)
colormap(gca,cmap)
%hcb2=colorbar;
pbaspect([1.333,1,1])

set(hsub(1:4),'XTickLabel',''); set(hsub,'YTickLabel','')

if printing==1
    cd(InputUser.ResultsDir)
    print(gcf,['WatershedCheck'],'-dpng','-r300');
end

clear hcb2 hcb1 hsub wsub cmap