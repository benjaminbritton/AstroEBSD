function FolderOut=Folder_Run(InputUser,Crystal_UCell,Crystal_LUT,...
    Settings_PCin,Settings_Cor,Settings_Rad,Settings_LUT,Phase_Num,pattern_list,ClockStart)

%run the folder
num_patterns=numel(pattern_list);

PC_out=zeros(3,Phase_Num,num_patterns);
PC_err=zeros(Phase_Num,num_patterns);
PC_errmax=zeros(num_patterns,1);
Phase=zeros(num_patterns,1);
PC_f=zeros(3,num_patterns);
rotdata=cell(num_patterns,Phase_Num);
banddata=cell(num_patterns,Phase_Num);
Peak_Centres=cell(num_patterns,1);

pTime('Radon Transform & Pattern Read',ClockStart);


parfor n=1:num_patterns
    %read the pattern
    [PatternIn] = ReadEBSDFile(pattern_list{n},InputUser.PatternFlip);
    
    %correct & convert
    [PatternCor] = EBSP_BGCor(PatternIn,Settings_Cor );
    % radon convert & Peak ID
    [ Peak_Centres{n}] = EBSP_RadHunt(PatternCor,Settings_Rad);
    
end

%do for 1st pattern for settings
[PatternIn] = ReadEBSDFile(pattern_list{1},InputUser.PatternFlip);
[~,PatternInfo ] = EBSP_BGCor(PatternIn,Settings_Cor );

pTime('Finding Pattern Centres',ClockStart);

parfor n=1:num_patterns
    PC_out_t=[];
    if InputUser.PCSearch == 0
        for num_P=1:Phase_Num
            PC_out_t(:,num_P)=Settings_PCin.start;
        end
        PC_out(:,:,n)=PC_out_t;
        PC_f(:,n)=PC_out_t;
    else
        
        %set up the GA
        PC_GA_options = optimoptions('ga');
        PC_GA_options.FunctionTolerance=1E-3;
        PC_GA_options.UseParallel=0;
        PC_GA_options.MaxGenerations=15;
        PC_GA_options.PopulationSize=30;
        PC_GA_options.MaxStallGenerations=20;
        PC_GA_options.Display='off';
        
        PC_GA_ub=Settings_PCin.start+Settings_PCin.range;
        PC_GA_lb=Settings_PCin.start-Settings_PCin.range;
        
        PC_err_t=zeros(Phase_Num,1);
        PC_out_t=zeros(3,Phase_Num);
        
        for num_P=1:Phase_Num
            FitFunc = @(PC_test) PC_GAOpt( PC_test,Peak_Centres{n},PatternInfo.size,Crystal_LUT,Crystal_UCell,num_P);
            [PC_out_t(:,num_P), PC_err_t(num_P)] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
        end
        
        PC_err(:,n)=PC_err_t;
        PC_out(:,:,n)=PC_out_t;
        
        [PC_errmax(n),Phase(n)]=nanmax(PC_err_t);
        PC_f(:,n)=PC_out_t(:,Phase(n));
    end
    
end

pTime('Indexing',ClockStart);

parfor n=1:num_patterns
    %convert the normals
    [nhat_gnom] = EBSP_NormConv(Peak_Centres{n},PatternInfo.size,PC_f(:,n)); %#ok<PFBNS>
    
    rotdata_one=cell(Phase_Num,1);
    banddata_one=cell(Phase_Num,1);
    
    %index for these phases
    for num_P=1:Phase_Num
        [rotdata_one{num_P},banddata_one{num_P}]=EBSP_Index(nhat_gnom,Crystal_LUT{Phase_Num},Settings_LUT{Phase_Num}.thresh_trig,Crystal_UCell{Phase_Num},eye(3)); %#ok<PFBNS>
    end
    
    rotdata(n,:)=rotdata_one;
    banddata(n,:)=banddata_one;
    
end

FolderOut.patternlist=pattern_list;
FolderOut.rotdata=rotdata;
FolderOut.banddata=banddata;
FolderOut.PC_errmax=PC_errmax;
FolderOut.PC_err=PC_err;
FolderOut.PC_out=PC_out;