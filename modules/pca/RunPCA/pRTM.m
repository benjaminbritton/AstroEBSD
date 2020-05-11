%% Plot results of RTM  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function []=pRTM(RTMOutput,tile,EDS,printing,InputUser,PCA_Setup)

InputPats=tile.pat;

figure('Position', get(groot,'ScreenSize'))
cmap1=cbrewer('qual','Set1',21);
cmap=repmat(cmap1,1333,1);
%cmap2=cbrewer('div','RdGy',100);
cmap3=cbrewer('div','PuOr',100);
cmap4=cbrewer('div','RdBu',100);
%figure('Position',[100,350,1500,600])
[hsub,wsub]=tight_subplot(2,3,0.03,0.03);

axes(hsub(1))
%subplot(2,3,1)
imagesc(RTMOutput.phase)
colormap(gca,cbrewer('qual','Dark2',10))
title('RTM phase assignment')
hcb1=colorbar;
pbaspect([1.33,1,1])
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
caxis(gca,[0,10]);

axes(hsub(2))
%subplot(2,3,2)
imagesc(RTMOutput.qual)
colormap(gca,'gray')
title('RTM quality map')
hcb2=colorbar;
pbaspect([1.33,1,1])
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
caxis(gca,[0,1]);

if exist(EDS)==1
axes(hsub(3))
%subplot(2,3,3)
imagesc(EDS.Maps(:,:,8))
colormap(gca,cmap3)
title('Mo peak height')
hcb2=colorbar;
pbaspect([1.33,1,1])
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
else
axes(hsub(3))
%subplot(2,3,3)
imagesc(RTMOutput.euler3)
colormap(gca,'parula')
title('Euler 1')
pbaspect([1.33,1,1])
hcb2=colorbar;
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
end

axes(hsub(4));
%subplot(2,3,4);
imagesc(tile.map_reshaped)
colormap(gca,cmap)
title('Principal components')
hcb2=colorbar;
pbaspect([1.33,1,1])
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')

if exist(EDS)==1
axes(hsub(5));
%subplot(2,3,5)
imagesc(EDS.Maps(:,:,3))
colormap(gca,cmap3)
title('Cr peak height')
pbaspect([1.33,1,1])
hcb2=colorbar;
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
else
axes(hsub(5))
imagesc(RTMOutput.euler2)
colormap(gca,'parula')
title('Euler 2')
pbaspect([1.33,1,1])
hcb2=colorbar;
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')    
end

axes(hsub(6));
%subplot(2,3,6)
imagesc(RTMOutput.euler1)
colormap(gca,'parula')
title('Euler 1')
pbaspect([1.33,1,1])
hcb2=colorbar;
set(gca,'XTickLabel',''); set(gca,'YTickLabel','')
%caxis(gca,[0,1]);

if printing==1
    cd(InputUser.ResultsDir)
    print(gcf,['RTM'],'-dpng','-r300');
end

clear hcb2 hcb1 hsub wsub cmap cmap1 cmap2 cmap3
end