%% Look for RCC no
subplot(1,2,1)
imagesc(tile.map_reshaped)
pbaspect([1.333,1,1])
subplot(1,2,2)
imagesc(RTM.Output.qual)
pbaspect([1.333,1,1])


%% Compare matched library patterns to RC-EBSPs
rows=4;
cols=4;
RCno=146;


imnos=length(InputUser.Phases)+1;

figure('Position', get(groot,'ScreenSize'))
[hsub,wsub]=tight_subplot(cols,rows,0.03,0.03);

for ind=1:imnos-1
    
    axes(hsub(ind))
    pat=reshape(RTM.matchedpats(ind,RCno,:,:),RTM.screensize,RTM.screensize);
    imagesc(pat)
    pbaspect([1,1,1])
    set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
    title(InputUser.Phases(ind))
end

axes(hsub(ind+1))
pat=reshape(tile.pat(:,:,RCno),RTM.screensize,RTM.screensize);
imagesc(pat)
pbaspect([1,1,1])
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
title('RCEBSP')

colormap('gray')

%%
figure
pat=reshape(tile.pat(:,:,RCno),RTM.screensize,RTM.screensize);
subplot(1,2,1)
a=imagesc(pat)
colormap(gca,'gray')
caxis(gca,[-15,15]);
pbaspect([1,1,1])

pat=reshape(RTM.matchedpats(2,RCno,:,:),RTM.screensize,RTM.screensize);
subplot(1,2,2)
b=imagesc(pat)
caxis(gca,[-3,3])
pbaspect([1,1,1])
colormap(gca,'hot')




