% TiledPCA  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Now perform PCA etc for each tile of the AOI
for tiling_no=1:(PCA_Setup.crop_factor.^2) %will =1 if run_one_tile enabled
disp(['Running tile #',num2str(tiling_no),' of ' num2str((PCA_Setup.crop_factor.^2))]);

%grab the correct tile if run_one_tile exists
if isnan(PCA_Setup.run_one_tile)~=1
    PCA_Setup.crop_factor=savetiling;
    tiling_no=PCA_Setup.run_one_tile;
    clear savetiling
end

%% Load the data
[Data,tile]=loadtile(tile,tiling_no,EBSD_Info,MapInfo,PCA_Setup,Settings_Cor,t1); %loads patterns and spectra and normalises to intensity and stdev independently

%% Set up the test matrix

%Go to the tile's corresponding folder
if printing==1
    cd(InputUser.ResultsDir)
    cd(['Tile',num2str(tiling_no)])
    printingfol=pwd;
end

%Grab the correct array for a tile in the AOI
if PCA_Setup.PCA_EBSD==1 & PCA_Setup.PCA_EDX==1
    testArray=[Data.Patterns_reshaped_norm;Data.EDSData_cor_normvector];
end

if PCA_Setup.PCA_EBSD==1 & PCA_Setup.PCA_EDX==0
    testArray=[Data.Patterns_reshaped_norm];
end

if PCA_Setup.PCA_EDX==1 & PCA_Setup.PCA_EBSD==0
    testArray=[Data.EDSData_cor_normvector];
end

%Grab the cropped patterns and SEM images
Plots.RadonCropped(:,:,tiling_no)=Plots.Radon(tile.rowstart(tiling_no):tile.rowfin(tiling_no), tile.colstart(tiling_no):tile.colfin(tiling_no));
Plots.SEMImageCropped(:,:,tiling_no)=Plots.SEMImage(tile.rowstart(tiling_no):tile.rowfin(tiling_no), tile.colstart(tiling_no):tile.colfin(tiling_no));

% normalise standard deviation for WHOLE testArray.
%testArray=testArray./std(testArray);
%remove nan values
testArray(isnan(testArray))=0;
% Clear patterns that aren't needed
clear Data 
    
%% Initialise 'Output' structure which stores PCA (unrotated) results.
Output.weighting_ratio=weighting;

% Run the PCA on the testArray found for this tile of the AOI
pTime('Running PCA',t1);
rPCA_aoi; %This is independent of what is contained in testArray

% Save a couple of bits of variance that are interesting...
Variance.explained=Output.explained;
Variance.latent=Output.latent;

% Print non-PCA figures for this AOI tile (image, IQ...)
if printing==1
    figure;
    imagesc(Plots.SEMImageCropped(:,:,tiling_no)); axis image off; colormap('gray'); title('SEM Image');
    print(gcf,[InputUser.name,'SEM_cropped'],'-dpng','-r300');

    figure
    imagesc(Plots.RadonCropped(:,:,tiling_no));axis image off;colormap('gray'); title('IQ');
    print(gcf,[InputUser.name,'_PQcropped'],'-dpng','-r300');
end

%% VARIMAX rotation
% select components using a watershed transform. (uncomment this if you'd
% like to see the result...)
%SelectComps_Watershed
%ret_comps=tile.Watershed_selcomps;

%Select components with variance tolerance
ret_comps=max(find(Output.explained>vartol));

if isnan(PCA_Setup.components)~=1;
    ret_comps=PCA_Setup.components;
end

%Needs to be at least two for VM rotation
if ret_comps<2;
    ret_comps=2;
end
tile.ret_comps(tiling_no)=ret_comps;

%% Reduce and rotate with optimised number of components
pTime(['Tile ',num2str(tiling_no),' selected ',num2str(ret_comps),' components.'],t1);
rPCA_aoi_VMrotation

VMOutput.explained=Output.explained;

if printing==1
    pPCA_aoi_VM;
end

%% Get averages
%Get average spectra
if PCA_Setup.PCA_EDX==1
    av_spectra=zeros(2048,ret_comps);
    for ind=1:ret_comps;
        
        %location in full map assigned to this component number
        [fullmap_loc]=find(VMOutput.PCA_VM_num==ind);
        
        %preallocate locations
        locs=zeros(size(fullmap_loc));
        
        %cycle through the locs and grab where that spatial location was in
        %the testarray 
        for l=1:length(fullmap_loc)
            locs(l)=find(fullmap_loc(l)==tile.xy_loc);
        end
        
        rawspectra=testArray(end-2047:end,locs);
        av_spectra(:,ind)=mean(rawspectra,2);
        clear rawspectra locs
    end
end

%assign (and concatenate) useful bits to the tile structure
if PCA_Setup.PCA_EDX==1
    tile.av_spectra=cat(2,tile.av_spectra,av_spectra);
    tile.spec=cat(2,tile.spec,VMOutput.spectrum_VM_n);
end

clear av_spectra

if PCA_Setup.PCA_EBSD==1
    tile.pat=cat(3,tile.pat,VMOutput.PCA_VM_Pat_n);
end

tile.map(:,:,tiling_no)=VMOutput.PCA_VM_num;
tile.filtmap(:,:,tiling_no)=VMOutput.filtered_map;
tile.dims(:,tiling_no)=[VMOutput.width,VMOutput.height];
tile.(['Tile',num2str(tiling_no)])=VMOutput;
tile.ret_comps(tiling_no)=VMOutput.NumPCA;
tile.variance=Variance;
tile.loadingmap(:,:,tiling_no)=VMOutput.PCA_VM_val;

close all
clear VMOutput Output testArray Variance printingfol
end
pTime(['Tiling complete.'],t1);
clear ret_comps tiling_no ind