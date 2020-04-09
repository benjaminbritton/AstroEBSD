%% Plot correlation similarities - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

figure('Position', get(groot,'ScreenSize'))
hold on
subplot(1,3,1)
EDSWindowing.qualmap=zeros(size(tile.map_reshaped));
for a=1:total_ret_comps
    locs=find(tile.map_reshaped==a);
    EDSWindowing.qualmap(locs)=EDSWindowing.EDScorrqual(a);
    clear locs
end
imagesc(EDSWindowing.qualmap)
pbaspect([1,1,1])
title('Quality');
cmap=cbrewer('div','PRGn',100);
colormap(cmap)
%
subplot(1,3,2)
EDSWindowing.skipmap=zeros(size(tile.map_reshaped));
for a=1:total_ret_comps
    locs=find(tile.map_reshaped==a);
    EDSWindowing.skipmap(locs)=EDSWindowing.skipped(a);
    clear locs
end
imagesc(EDSWindowing.skipmap)
pbaspect([1,1,1])
title('Number skipped');
cmap=cbrewer('div','RdBu',100);
colormap(cmap)

subplot(1,3,3)
EDSWindowing.countmap=zeros(size(tile.map_reshaped));
for a=1:total_ret_comps
    locs=find(tile.map_reshaped==a);
    EDSWindowing.countmap(locs)=EDSWindowing.count(a);
    clear locs
end
imagesc(EDSWindowing.countmap)
pbaspect([1,1,1])
title('Dissimilar number counted');
cmap=cbrewer('div','RdBu',100);
colormap(cmap)
%colorbar
hold off

if printing==1
    cd(InputUser.ResultsDir)
    print(gcf,['EDS_Correlation'],'-dpng','-r300');
end


clear cmap 