% pPCA_aoi_VM  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% This plots the PCA datasets and prints them as bmps

cd(printingfol);
ItFolder=strrep(['RotPCA','_n',num2str(VMOutput.NumPCA),'_w',num2str(Output.weighting_ratio,'%.2e')],'.',',');
mkdir(ItFolder);
cd(ItFolder);

%% Plot the PCs
%plot the new factors
figure('Position', [400, 300, 900, 600]);
for n=1:VMOutput.NumPCA
        
if PCA_Setup.PCA_EBSD==1;
    subplot('Position',[0,0.25,0.5,0.5]);
    imagesc(VMOutput.PCA_VM_Pat_n(:,:,n)); axis image; axis xy; axis tight; colormap('gray'); axis off; title(['Kikuchi ',int2str(n)]);
end
        
if PCA_Setup.PCA_EDX==1
    subplot('Position',[0.55,0.35,0.4,0.35]);
    plot(VMOutput.spectrum_VM_n(:,n)); title(['Spectrum ', int2str(n)]);
    xlim([0, 1200]);
    if max(VMOutput.spectrum_VM_n(:,n))==0
    else   
      ylim([0, 1.01*max(VMOutput.spectrum_VM_n(:,n))]); %max y value normalised to this spectrum
    end 
      %very crude EDX fingerprint...
%     subplot('Position',[0.55,0.1,0.4,0.3]);
%     [pks,locs]=findpeaks(VMOutput.spectrum_VM_n(:,n),'MinPeakProminence',0.0009);
%     pks_norm=pks./sum(pks);
%     plot(pks_norm); xlim([0,15]); ylim([0,0.25]); title(['EDX fingerprint ', int2str(n)]);
%     xlabel('Peak number'); ylabel('Normalised peak intensity');
%     
    if printing==1;
        print(gcf,[InputUser.name,'_rPCAfactors_',num2str(n)],'-dpng','-r300');
    end
end
end
clear pks locs pks_norm n

figure
for n=1:VMOutput.NumPCA
    %subplot(PCAGrid(1),PCAGrid(2),n);
    imagesc(VMOutput.PCA_VM_Map_n(:,:,n)); axis image; axis ij; axis tight; axis off;title(int2str(n));
    
    if printing==1
    print(gcf,[InputUser.name,'_rPCAfactorsmap_',num2str(n)],'-dpng','-r300');
    end
end
clear n 

%plot the factor map
figure;
imagesc(VMOutput.PCA_VM_num); axis image; colormap('colorcube');
if printing==1
print(gcf,[InputUser.name,'_rPCA'],'-dpng','-r300');
end

%plot the quality map
figure;
imagesc(VMOutput.PCA_VM_val); axis image; colormap('gray');
if printing==1
print(gcf,[InputUser.name,'_PCAq'],'-dpng','-r300');
end

% scree plot: components vs proportion of total variance retained by each
figure
plot(Output.explained);
xlim([1,VMOutput.NumPCA]);
xlabel('Principal component');
ylabel('Percentage of variance retained');
if printing==1
print(gcf,[InputUser.name,'Scree_',num2str(VMOutput.NumPCA)],'-dpng','-r300');
end


if PCA_Setup.weighted==1
    save([InputUser.name,'_PCA_','w',num2str(Output.weighting_ratio),'_n',num2str(VMOutput.NumPCA),'_wkspc',],'VMOutput','InputUser','MicroscopeData','PCA_Setup','Plots','Settings_Cor')
else 
    save([InputUser.name,'_PCA_','_wkspc',],'VMOutput','InputUser','MicroscopeData','PCA_Setup','Plots','Settings_Cor')
end
cd(InputUser.ResultsDir)

close all

clear ItFolder





