%Readout settings and check for incompatabilities in settings

disp(' ')
disp(' ')
disp('* * * * * * * * * * * * * * * * * * * * * * * * * * * *')

%% Check for incompatibilities in settings

% run_one_tile > crop_factor
if PCA_Setup.run_one_tile>PCA_Setup.crop_factor^2
    error("Can't run a single tile with position more than crop_factor squared! \n Edit PCA_Setup.run_one_tile to be less than (PCA_Setup.crop_factor)^2",class(1))
end

% kernel size not odd
if PCA_Setup.SpatialKernel==1 & iseven(PCA_Setup.KernelRadius)
    error("Kernel radius must be odd")
end

if PCA_Setup.PCA_EBSD==0 & PCA_Setup.weighted~=0
    PCA_Setup.weighted=0;
    warning('Not running EBSD, so PCA_Setup.weighted changed to zero')
end

if PCA_Setup.PCA_EDX==0 & PCA_Setup.weighted~=0
    PCA_Setup.weighted=0;
    warning('Not running EDS, so PCA_Setup.weighted changed to zero')
end

binaryvalues=[PCA_Setup.PCA_EBSD,PCA_Setup.PCA_EDX,PCA_Setup.weighted,PCA_Setup.RTM,PCA_Setup.MedFilt,PCA_Setup.SpatialKernel,printing,Settings_Cor.gfilt,Settings_Cor.radius,Settings_Cor.resize,Settings_Cor.SplitBG,Settings_Cor.RealBG,Settings_Cor.Square,Settings_Cor.SquareCrop,Settings_Cor.LineError,Settings_Cor.MeanCentre];
binaryargs={'PCA_Setup.PCA_EBSD','PCA_Setup.PCA_EDX','PCA_Setup.weighted','PCA_Setup.RTM','PCA_Setup.MedFilt','PCA_Setup.SpatialKernel','printing','Settings_Cor.gfilt','Settings_Cor.radius','Settings_Cor.resize','Settings_Cor.SplitBG','Settings_Cor.RealBG','Settings_Cor.Square','Settings_Cor.SquareCrop','Settings_Cor.LineError','Settings_Cor.MeanCentre'};

for i = 1:length(binaryvalues)
    if binaryvalues(i)~=1 & binaryvalues(i)~=0
        error(['Variable ',binaryargs{i},' should be binary'])
    end
end

%% Tell the user what they're running

if PCA_Setup.RTM==1
    RTM_str=' with RTM ';
else
    RTM_str=' ';
end

if PCA_Setup.PCA_EBSD==1
    data_str=' EBSD';
    if PCA_Setup.PCA_EDX==1;
        data_str=[data_str,' and EDS'];
    end
else
    data_str=' EDS ';
end

if PCA_Setup.crop_factor>1 & isempty(PCA_Setup.run_one_tile)==1
    tiling_str=[', with ',num2str(PCA_Setup.crop_factor^2),' tiles.'];
    warn_tiles=0;
    
elseif PCA_Setup.crop_factor>1 & isempty(PCA_Setup.run_one_tile)==0
    tiling_str=[', on a single tile (',num2str(PCA_Setup.run_one_tile),' of ',num2str(PCA_Setup.crop_factor^2),')'];
    warn_tiles=1;
    
else
    tiling_str=[];
    warn_tiles=2;
end

DataTypeStr=['Running PCA',RTM_str,'on',data_str,' data',tiling_str];
disp(DataTypeStr)

if warn_tiles==1
    warning('You are running on one tile of this dataset')
end

if warn_tiles==2
    warning('You are not using any tiling. Be careful of RAM issues.')
end

%% What parameters?

if PCA_Setup.PCA_EBSD==1 & PCA_Setup.PCA_EDX==1 & PCA_Setup.weighted==1
    weighting_str=[' EBSD weighting ',num2str(PCA_Setup.EBSD_weighting),' and'];
    warn_weighting=0;
elseif PCA_Setup.PCA_EBSD==1 & PCA_Setup.PCA_EDX==1 & PCA_Setup.weighted==0
    weighting_str=[' no EBSD weighting and'];
    warn_weighting=1;
else
    weighting_str=[];
    warn_weighting=0;
end

vartol_str=num2str(PCA_Setup.variance_tolerance);

SettingsStr=['...using',weighting_str, ' variance tolerance ',vartol_str];

disp(SettingsStr)

if warn_weighting==1
    warning('Are you sure you want to run this correlative dataset with no weighting?')
end

if PCA_Setup.SpatialKernel==1
    KernelStr=['...and a spatial kernel of size ',num2str(PCA_Setup.KernelRadius)];
    disp(KernelStr)
end

if PCA_Setup.RTM==1
    disp(['...with template matching for ',num2str(length(InputUser.Phases)),' phases, at SO3 sampling frequency of ',num2str(RTM_setup.Sampling_Freq) ' degrees'])
end

if isempty(PCA_Setup.components)==0
    warning(['You are running with ',num2str(PCA_Setup.components),' components fixed for every tile']);
end

disp('* * * * * * * * * * * * * * * * * * * * * * * * * * * *')
disp(' ')
disp(' ')
