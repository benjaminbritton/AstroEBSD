function [EBSP_Fe_PCOut] = EBSP_PCSearch(EBSP_Fe_PC,Settings_Cor,Settings_Rad,Settings_PCin,Phase_Num,Crystal_LUT,Crystal_UCell)
%EBSP_PCSEARCH Find the pattern centre

%set the bounds on the GA
PC_GA_ub=Settings_PCin.start+Settings_PCin.range;
PC_GA_lb=Settings_PCin.start-Settings_PCin.range;


%set up the GA - needs the global optimisation toolbox
PC_GA_options = optimoptions('ga');
PC_GA_options.FunctionTolerance=1E-3;
PC_GA_options.UseParallel=0;
PC_GA_options.MaxGenerations=15;
PC_GA_options.PopulationSize=30;
PC_GA_options.MaxStallGenerations=20;
PC_GA_options.Display='iter';
%for help with this:
%doc ga

EBSP_Fe_PC.PC_out=zeros(3,Phase_Num);
EBSP_Fe_PC.PC_err=zeros(Phase_Num,1);

%Radon transform
[ EBSP_Fe_PC.Peak_Centre,EBSP_Fe_PC.Single.Peak_Set_All,EBSP_Fe_PC.Peak_Set_All,...
            EBSP_Fe_PC.R_EBSP,EBSP_Fe_PC.R_Edge,EBSP_Fe_PC.R_rho,EBSP_Fe_PC.R_theta ] ...
            = EBSP_RadHunt( EBSP_Fe_PC.PatternIn,Settings_Rad);
        
%index & PC hunt for all the phases
for num_P=1:Phase_Num
    FitFunc = @(PC_test) PC_GAOpt( PC_test,EBSP_Fe_PC.Peak_Centre,Settings_Cor.size,Crystal_LUT,Crystal_UCell,num_P);
    [EBSP_Fe_PC.PC_out(:,num_P), EBSP_One.PC_err(num_P)] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
end
%find the best fitting
[EBSP_Fe_PC.PC_errmax,EBSP_Fe_PC.Phase]=nanmax(EBSP_Fe_PC.PC_err);

%assign that PC
EBSP_Fe_PC.PC=EBSP_Fe_PC.PC_out(:,EBSP_Fe_PC.Phase);

%provide the outputs
EBSP_Fe_PCOut=EBSP_Fe_PC;

end

