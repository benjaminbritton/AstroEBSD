%% Tile reassembly  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Reassemble from tile AOIs into full maps
if isnan(PCA_Setup.run_one_tile)~=1
    tile.map_reshaped=tile.map(:,:,PCA_Setup.run_one_tile);

    if PCA_Setup.PCA_EDX==1
        tile.spec_reshaped=reshape(tile.spec,2048,sum(tile.ret_comps(1:end)));
    end

    if PCA_Setup.PCA_EBSD==1
        tile.pat_reshaped=reshape(tile.pat,EBSD_Info.PatSizeH*EBSD_Info.PatSizeW,sum(tile.ret_comps(1:end)));
    end
else
    CombineTiles
end

%Median filter the outputted map (can be handy... don't have to use)
tile.map_reshaped_original=tile.map_reshaped;

if PCA_Setup.MedFilt==1
    tile.map_reshaped=medfilt2(tile.map_reshaped);
end
    
%Print assignment map
if printing==1
    cd(InputUser.ResultsDir)
    cmap1=cbrewer('qual','Set1',80);
    cmap=repmat(cmap1,100,1);
    figure
    imagesc(tile.map_reshaped); colormap(cmap); axis off;
    print(gcf,[InputUser.name,'_AssignmentMap_','tile_',num2str(PCA_Setup.crop_factor^2)],'-dpng','-r300');
    clear cmap1 cmap
end

%Plot correlations of EDS spectra
if PCA_Setup.PCA_EDX==1
    %EDS_PC_correlationcount;
    %EDS_PC_correlationcount_plotting;
    %ElementMap; %This was for a specific use case - may be interesting to
    %adapt.
end
