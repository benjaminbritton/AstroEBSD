%% The CODE  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%cd to the HDF5 folder
cd(InputUser.HDF5_folder)
PCA_Setup.cropping=1;


% initialise orientation matrices
RTM_setup.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM_setup.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM_setup.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

% Turn off warnings for breadHDF5
InputUser.breadhdf5_warnings_off=1;

%% Sets up required plugins, loads data and calculates cropping windows.
for fileloop=1:length(files)

clear tile RTM Plots EDSWindowing
RTM=RTM_setup;

InputUser.HDF5_file = files{fileloop};  %h5 file, include the file extension
t1=clock;
if exist(InputUser.SavingFolder) == 0
    mkdir(InputUser.SavingFolder);
end

cd(InputUser.SavingFolder);

%% Creates appropriate filename and HDF5 handling parameters.
PCA_filenamesetup;

%% Load metadata
[MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadHDF5( InputUser );

[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);

Plots.Radon=MapInfo.Data_InputMap.RadonQuality;
Plots.SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));

MapInfo.MapSize1=max(MapInfo.MapData.YBeam)-min(MapInfo.MapData.YBeam)+1; % number of rows
MapInfo.MapSize2=max(MapInfo.MapData.XBeam)-min(MapInfo.MapData.XBeam)+1; % number of columns

%% Refine the pattern centre if requested
if Refine.run==1
    [MapInfo]=fRefineMap(Refine,RTM,MapInfo,MicroscopeData,InputUser,Settings_Cor,t1);
end

%% Plot an example pattern
figure;
[ RefPat ] = bReadEBSP(MapInfo.EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
subplot(1,2,1); imagesc(RefPat); axis image; axis xy; axis tight; colormap('gray'); title('Input pattern');
subplot(1,2,2); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');

EBSD_Info.PatSizeH=size(RefPat_cor,2);
EBSD_Info.PatSizeW=size(RefPat_cor,1);
clear RefPat RefPat_cor

%% Get tiling parameters
MapInfo=TilingSetup(MapInfo,PCA_Setup.crop_factor);

%% Now loop over the inputted weights

for weightloop=1:length(PCA_Setup.EBSD_weighting)%% Setup PCA, load data, generate tiles.

%Main loop over tiled areas
%Results from full AOI are stored in 'tile' structure.
%Results from RTM are stored in 'RTM' structure.
clear tile

vartol=PCA_Setup.variance_tolerance(1);
weighting=PCA_Setup.EBSD_weighting(weightloop);
disp(['weighting  ',num2str(weightloop),' of ',num2str(length(PCA_Setup.EBSD_weighting))]);
cd(InputUser.SavingFolder)

%pre-allocate empty arrays for spectra, patterns, etc for concatenation.
tile.pat=[];
tile.spec=[];
tile.av_spectra=[];
tile.weight=PCA_Setup.EBSD_weighting(weightloop);

if printing==1
    FolderSetup_PCA;
end
close all

%do a quick adjustment of cropfactor if only want to run one tile
if isnan(PCA_Setup.run_one_tile)~=1
    savetiling=PCA_Setup.crop_factor;
    PCA_Setup.crop_factor=1;
end

%% Analyse the data

TiledPCA;
TiledPCA_Reassembly;

%Run the RTM if required
if PCA_Setup.RTM==1
   RTM_PCA;
   %RTM.Output=RTMOutput;
   clear RTMOutput;
end

%% Save results
%Save the results. Will overwrite as it goes through the
%weightloops/tolerances...
if printing==1
    cd(InputUser.ResultsDir)
    if PCA_Setup.RTM==1
        if PCA_Setup.PCA_EDX==1
            save(['Results'],'PCA_Setup','tile','RTM','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        else
            save(['Results'],'PCA_Setup','tile','RTM','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        end
    
    else 
        if PCA_Setup.PCA_EDX==1
            save(['Results'],'PCA_Setup','tile','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        else
            save(['Results'],'PCA_Setup','tile','InputUser','weighting','Settings_Cor','MicroscopeData','EBSD_Info','MapInfo','-v7.3')
        end
    end   
end

%% Spectra saving for RC-spectra
if PCA_Setup.PCA_EDX==1
    if printing==1
        cd(InputUser.ResultsDir)
        mkdir('RCSpectra')
        cd('RCSpectra')

        %Grabs a .spx file from somewhere else, then changes channel numbers to
        %those of RC-spectra. NB operates on tile.spec if needs to change
        SPX_replace_PCA(fullfile(InputUser.Astro_loc,'modules','pca','Analysis','Bruker_Spectra_Hack'),'spec.spx',tile.spec,pwd)
        cd(InputUser.ResultsDir)

        mkdir('AverageSpectra')
        cd('AverageSpectra')
        %Grabs a .spx file from somewhere else, then changes channel numbers to
        %those of RC-spectra. NB operates on tile.spec if needs to change
        SPX_replace_PCA(fullfile(InputUser.Astro_loc,'modules','pca','Analysis','Bruker_Spectra_Hack'),'spec.spx',tile.av_spectra,pwd)
        cd(InputUser.ResultsDir)

    end

    clear ret_comps sel_comps weighting ComponentMap ElementMaps InputPats tiling_no total_ret_comps n
    %clear tile EDSWindowing
end

clear weighting vartol
end %end of weightloop
clear weightloop

cd(InputUser.SavingFolder)
close all

pTime(['completed ',num2str(fileloop), ' of ' num2str(length(files))],t1);
end %end of fileloop
clear fileloop ans
