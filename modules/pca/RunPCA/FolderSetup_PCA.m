% FolderSetup_PCA  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Generates folder name and prints SEM, Radon, EDS integrations for each
% tile.

% Generate folder name

cd(InputUser.SavingFolder)
vtol=PCA_Setup.variance_tolerance;

HomeDir=[InputUser.name,'_',PCA_Setup.suff];

if PCA_Setup.weighted==1 & PCA_Setup.SpatialKernel==1
    ResultsDir=['weighting',strrep(num2str(weighting),'.',','),'_vt',strrep(num2str(vtol),'.',','),'_',PCA_Setup.KernelFnLabel,'_rad',num2str(PCA_Setup.KernelRadius)];
elseif PCA_Setup.weighted==1 & PCA_Setup.SpatialKernel==0
    ResultsDir=['weighting',strrep(num2str(weighting),'.',','),'_vt',strrep(num2str(vtol),'.',',')];
elseif PCA_Setup.weighted==0 & PCA_Setup.SpatialKernel==1
    ResultsDir=['vt',strrep(num2str(vtol),'.',','),'_',PCA_Setup.KernelFnLabel,'_rad',num2str(PCA_Setup.KernelRadius)];
elseif PCA_Setup.weighted==0 & PCA_Setup.SpatialKernel==0
    ResultsDir=['vt',strrep(num2str(vtol),'.',',')];
end
    
clear suff vtol
mkdir(HomeDir);
cd(HomeDir);
InputUser.HomeDir=pwd;

mkdir(ResultsDir)
cd(ResultsDir)
InputUser.ResultsDir=pwd;

clear ResultsDir


% Print SEM, Radon and EDS integration for full map
figure;
imagesc(Plots.SEMImage); axis image off; colormap('gray'); title('SEM Image');
print(gcf,[InputUser.name,'SEM'],'-dpng','-r300');

figure
imagesc(Plots.Radon);axis image off;colormap('gray'); title('IQ');
print(gcf,[InputUser.name,'_PQ'],'-dpng','-r300');

% if PCA_Setup.PCA_EDX==1;
%    figure
%    imagesc(Plots.EDSInteg)
%    print(gcf,[name,'EDS_integ'],'-dpng','-r300'); 
% end

%% Make all tile folders
for i=1:(PCA_Setup.crop_factor.^2);
    mkdir(['Tile',num2str(i)]);
    cd(['Tile',num2str(i)]);
    cd(InputUser.ResultsDir);
end

clear i HomeDir