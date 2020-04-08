%% Pat comparison script  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% if RTM_update.updatelocs==1
%     %fix for having changed phase folders
%     RTM.Bin_loc=RTM_update.Bin_loc;
%     RTM.Phase_Folder=RTM_update.Phase_folder;
% end

%preallocate zeros
matchedpat=zeros(EBSD_Info.PatSizeH,EBSD_Info.PatSizeW);


figure('Position', get(groot,'ScreenSize'))
subplot(subplot_h,subplot_w,1)
imagesc(InputPats(:,:,component))
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('PCA pattern')
pbaspect([1,1,1])

for iii=1:length(InputUser.Phases)
subplot(subplot_h,subplot_w,iii+1)

[RefinedTemplatePat_uncor]=RefineTemplate(iii,component,RTM,tile,MapInfo,MicroscopeData,EBSD_Info,Settings_Cor,InputUser);

RefinedTemplatePat=EBSP_BGCor(RefinedTemplatePat_uncor,Settings_Cor);

imagesc(RefinedTemplatePat)
title([num2str(iii), InputUser.Phases(iii),num2str(RTM.Output.PeakHeight(component,iii))])
pbaspect([1,1,1])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
end
colormap('gray')

%%
if patcomp.print==1

    % Print comparison patterns to current directory
    %clear matchedpat component iii
    %title(['Phase: ', InputUser.Phases(phase), ' / Quality: ',num2str(RTM.Output.peakheight(phase,component))])
    figure
    imagesc(InputPats(:,:,component))
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    pbaspect([1,1,1])
    colormap('gray')
    print(gcf,'pcapat','-dpng','-r300')

    for iii=1:length(InputUser.Phases)
    
    [RefinedTemplatePat]=RefineTemplate(iii,component,RTM,tile,MapInfo,MicroscopeData,EBSD_Info,Settings_Cor,InputUser);
    
    imagesc(RefinedTemplatePat)
    title([num2str(iii), InputUser.Phases(iii),num2str(RTM.Output.PeakHeight(iii,component))])
    pbaspect([1,1,1])
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap('gray')

    print(gcf,['pat',num2str(iii)],'-dpng','-r300')

    end
    close all
end

clear iii