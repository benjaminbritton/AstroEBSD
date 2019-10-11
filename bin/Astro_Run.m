%Pattern Indexing - Map Form
%This code should:
%(1) File handle the input file --> converting to suitable HDF5 if needed
%(2) Radon transform patterns in the map
%(3) Pattern centre search
%(4) Index
%(5) Data plots

try
    
    %% start the code
    %start the clock
    ClockStart=clock;
    
    %% determine the mode
    
    %InputUser.Mode = 'Isolated'; %modes =
    %Isolated = single file
    %Map_Single = single from map
    %Map_All = full map, with area PC search
    %Folder = blind folder
    
    
    if strcmpi(InputUser.Mode,'Isolated')
        Settings_Mode = 1;
    elseif strcmpi(InputUser.Mode,'Folder')
        Settings_Mode = 2;
    elseif strcmpi(InputUser.Mode,'Map_Single')
        Settings_Mode = 3;
    elseif strcmpi(InputUser.Mode,'Map_All')
        Settings_Mode = 4;
    else
        error(['The input setting mode ' InputUser.Mode ' is not supported']);
    end
    %% Build the phase data
    pTime('Building Phases',ClockStart);
    [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );
   
    

    %% BCF/Map loading
    
    if Settings_Mode == 3 || Settings_Mode == 4 % 3 Map_Single, 4 = Map_All
        pTime('BCF/HDF5 Loading',ClockStart);
        % convert the file to HDF5 & read data
        %read & build the HDF5
        [MapData,MicroscopeData,Phase_BrukerData,EBSD_DataInfo] = BCF_HDF5( InputUser );
        
        %read the map & convert to area data
        [Data_InputMap] = EBSD_Map(MapData,MicroscopeData);
        
        % generate the static BG
        pTime('Static background generation (if needed)',ClockStart);
        [ Settings_Cor ] = EBSP_StaticBG( Settings_Cor,MicroscopeData,EBSD_DataInfo);
    end
    
    %% Single pattern function
    if Settings_Mode == 1 %single pattern
        %load the pattern
        EBSP_One.PatternIn=ReadEBSDFile(InputUser.PatternLoc,InputUser.PatternFlip);
    end
    
    if Settings_Mode == 3 || Settings_Mode == 4  %load from map
        %load the pattern
        EBSP_One.P_num=Data_InputMap.PMap(InputUser.OnePatternPosition(2),InputUser.OnePatternPosition(1)); 
        EBSP_One.PatternIn=bReadEBSP(EBSD_DataInfo,EBSP_One.P_num);
    end
    
    %set up folder
        if Settings_Mode ==  2
        
        pTime('Folder Listing',ClockStart);
        [pattern_list,num_patterns]=Folder_Prep(InputUser);
        n=1;
        [EBSP_One.PatternIn] = ReadEBSDFile(pattern_list{n},InputUser.PatternFlip);
        
        end
    
    %adjust the pattern settings in GUI 
    
    %note that the BG correction for multi patterns is not working fully -
    %TBB 24/04/2017
    pTime('Review Settings (interactive - please close window to continue)',ClockStart);
    Astro_EBSPset(EBSP_One,Settings_Cor,Settings_Rad,Settings_PCin,InputUser)
    
    %% Full map functions
    
    if Settings_Mode == 4
        % Radon Transform for map
        pTime('Starting Radon Transforms',ClockStart);
        [ Peak_Centres,Peak_Quality,Peak_NBands,EBSD_Info ] = Map_Radon( Data_InputMap,EBSD_DataInfo,Settings_Cor,Settings_Rad );
        
        % Pattern centre search
        pTime('Starting Pattern Centre Search',ClockStart);
        [PCOut] = ...
            Map_PCSearch(Data_InputMap,Peak_Centres,Peak_NBands,EBSD_Info,Crystal_LUT,Crystal_UCell,Settings_PCin);
        
        % Index all patterns
        pTime('Starting Indexing Step',ClockStart);
        [Indexed_Rotdata,Indexed_Banddata] = ...
            Map_Index( Data_InputMap,Peak_Centres,Peak_NBands,Phase_Num,PCOut.Fit_2nd,Crystal_UCell,Crystal_LUT,Settings_LUT,EBSD_Info,MicroscopeData);
        
        % Convert Hough + rotation data into map form
        pTime('Map Generation',ClockStart);
        [Data_OutputMap] = Map_Generate( Data_InputMap,Indexed_Rotdata,Peak_Quality,Crystal_UCell,Phase_Num );
        
        %generate the quaternion maps
        Data_QMap=Map_Quats(Data_OutputMap);
        
        % Calculate IPF maps
        pTime('Generate IPF Maps',ClockStart);
        Plot_IPF.X = Map_IPF([1 0 0],Data_QMap);
        Plot_IPF.Y = Map_IPF([0 1 0],Data_QMap);
        Plot_IPF.Z = Map_IPF([0 0 1],Data_QMap);
    end
    
    %% Run a single pattern through
    if ismember(Settings_Mode,[1 3 4])
        
        [ EBSP_One.PatternCor,EBSP_One.PatternInfo ] = EBSP_BGCor( EBSP_One.PatternIn,Settings_Cor );
        % radon convert & Peak ID
        [ EBSP_One.Peak_Centre,EBSP_One.Single.Peak_Set_All,EBSP_One.Peak_Set_All,...
            EBSP_One.R_EBSP,EBSP_One.R_Edge,EBSP_One.R_rho,EBSP_One.R_theta ] ...
            = EBSP_RadHunt( EBSP_One.PatternCor,Settings_Rad);
        
        if InputUser.PCSearch == 1 && ismember(Settings_Mode,[1 2 3])
            
            %set up the GA
            PC_GA_options = optimoptions('ga');
            PC_GA_options.FunctionTolerance=1E-3;
            PC_GA_options.UseParallel=0;
            PC_GA_options.MaxGenerations=15;
            PC_GA_options.PopulationSize=30;
            PC_GA_options.MaxStallGenerations=20;
            PC_GA_options.Display='iter';
            
            PC_GA_ub=Settings_PCin.start+Settings_PCin.range;
            PC_GA_lb=Settings_PCin.start-Settings_PCin.range;
            
            EBSP_One.PC_out=zeros(3,Phase_Num);
            EBSP_One.PC_err=zeros(Phase_Num,1);
            
            for num_P=1:Phase_Num
                pTime(['Phase ' num2str(num_P) ' PC Search'],ClockStart);
                FitFunc = @(PC_test) PC_GAOpt( PC_test,EBSP_One.Peak_Centre,EBSP_One.PatternInfo.size,Crystal_LUT,Crystal_UCell,num_P);
                [EBSP_One.PC_out(:,num_P), EBSP_One.PC_err(num_P)] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
            end
            [EBSP_One.PC_errmax,EBSP_One.Phase]=nanmax(EBSP_One.PC_err);
                       
            EBSP_One.PC=EBSP_One.PC_out(:,EBSP_One.Phase);
        elseif InputUser.PCSearch == 0 && ismember(Settings_Mode,[1 2 3])
             EBSP_One.PC=Settings_PCin.start;
             
        elseif InputUser.PCSearch == 1 && Settings_Mode == 4
            EBSP_One.PC(1)=PCOut.Fit_2nd.PCx_map(InputUser.OnePatternPosition(2),InputUser.OnePatternPosition(1));
            EBSP_One.PC(2)=PCOut.Fit_2nd.PCy_map(InputUser.OnePatternPosition(2),InputUser.OnePatternPosition(1));
            EBSP_One.PC(3)=PCOut.Fit_2nd.PCz_map(InputUser.OnePatternPosition(2),InputUser.OnePatternPosition(1));
            
        end
        
        if isfield(EBSP_One,'PC') == 0
            EBSP_One.PC = Settings_PCin.start;
        end
        
            % Convert the bands to normal space
            [ EBSP_One.nhat_gnom] = EBSP_NormConv( EBSP_One.Peak_Centre,[EBSP_One.PatternInfo.size],EBSP_One.PC);
            
            %index for these phases
            for num_P=1:Phase_Num
            [EBSP_One.rotdata{num_P},EBSP_One.banddata{num_P}]=EBSP_Index(EBSP_One.nhat_gnom,Crystal_LUT{Phase_Num},Settings_LUT{Phase_Num}.thresh_trig,Crystal_UCell{Phase_Num},eye(3)); %#ok<PFBNS>
            end
            
            %generate the geometry
            [ EBSP_One.PatternGeometry ] = EBSP_Gnom( EBSP_One.PatternInfo,EBSP_One.PC );
    end
    

    
    if Settings_Mode == 2
        %run the folder
%         FolderOut=Folder_Run(InputUser,Crystal_UCell,Crystal_LUT,...
%             Settings_PCin,Settings_Cor,Settings_Rad,Settings_LUT,Phase_Num,pattern_list,pattern_list,ClockStart);
        FolderOut=Folder_Run(InputUser,Crystal_UCell,Crystal_LUT,...
            Settings_PCin,Settings_Cor,Settings_Rad,Settings_LUT,Phase_Num,pattern_list,ClockStart);

        %prepare one EBSP for plotting
        if InputUser.OnePatternPosition(1) > size(FolderOut.patternlist,1)
            InputUser.OnePatternPosition(1)=size(FolderOut.patternlist,1);
            disp('Pattern number for single plot reset to last pattern indexed');
        end
        
        EBSP_One=Folder_EBSPPrep(FolderOut,InputUser,Settings_Cor,Settings_Rad,InputUser.OnePatternPosition(1),1);
    end
    %%
    %save this data with a timestamp
    save(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_' pNameClock '.mat']));
    
    
    E_complete_text=pTime('AstroEBSD Index Complete',ClockStart);
 
catch err %% catch any errors tha happen
    
    % write the error to string for email
    % first line: message
    failed_time=pTime('Code Failed',ClockStart);
    E_message{1}=failed_time;
    E_message{2}=sprintf('%s\n',err.message);
    
    % following lines: stack
    for e=1:length(err.stack)
        E_message{e+2}=sprintf('%s in %s at %i\n',err.stack(e).name,err.stack(e).file,err.stack(e).line);
    end
    
  
    %save this data with a timestamp
    save(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_' pNameClock '.mat']));
        
    rethrow(err)
end
