function [ PCOut ] = ...
    Map_PCSearch(MapArea,Peak_Centres,Peak_Bands,outputs,Crystal_LUT,Crystal_UCell,PCin)
%PC_MAPSEARCH Search within a map for the pattern centre at each point
%(c) TBB 2017 - b.britton@imperial.ac.uk
% http://www.expmicromech.com
%
% Pattern Centre Searching Algorithm:
% uses two passes of a genetic algorithm for the search
% the genetic algorithms uses a weighted solution to the crystal orientation
% weighting is per the number of bands solved vs the total number input
% needs a seed pattern centre and a range of points to search
%
% Pattern centre data is output as per the conventions described in:
%
% T.B. Britton, J. Jiang, Y. Guo, A. Vilalta-Clemente, D. Wallis,
% L. Hansen, A. Winkelmann, A.J. Wilkinson 
% Tutorial: Crystal orientations and EBSD — Or which way is up? 
% Materials Characterization (2016)
% http://dx.doi.org/10.1016/j.matchar.2016.04.008
%
%
% INPUTS
%
%   EBSPData = EBSD map information for HDF5 reading
%   MapArea = Information on the Map
%   Settings_Cor = pattern background correction
%   Settings_Rad = radon transform settings
%   Crystal_LUT = structure containing look up table for phases
%   Crystal_UCell = unit cell information
%   PCin.array = [xpts, ypts] subset from the map
%   PCin.start = initial pattern centre guess - mean for map
%   PCin.range = range for 1st search, +/- value
%
% OUTPUTS
%
%   PC_data = Pattern Centre Data - 1st pass (subset)
%   PCFit_out = Fitted Data - 1st Pass (full map)
%   PC_data2 = Pattern Centre Data - 2nd pass (subset)
%   PCFit_out2 = Fitted Data - 2nd Pass (subset)


PC_subsize=50/1000; %pixel precision for second search, i.e. +/-0.05 pat frac

num_phases=size(Crystal_UCell,2);

%build the pattern centre grid
[PC_patnum,PC_patx,PC_paty,PC_num,PC_gridx,PC_gridy]=PC_GridBuild(MapArea,PCin.array);

%options - tuned as a balence between speed & accuracy 
%suitable for a 1st pass calculation
%with a follow up second calculation
%the population size and max generations are interesting variables
%population size = number of tries per generation - widen to deal with a
%rougher miniminsation
%max generations - number of generations allowed 
%i.e. how far does the search go
%
%note that the second stage of this will fit a grid, robustly, and use that
%to improve the results

PC_GA_options = optimoptions('ga');
PC_GA_options.FunctionTolerance=1E-3;
PC_GA_options.UseParallel=0;
PC_GA_options.MaxGenerations=15;
PC_GA_options.PopulationSize=30;
PC_GA_options.MaxStallGenerations=20;
PC_GA_options.Display='none';

%change the population size for pass 2
PC_GA_options2=PC_GA_options;
PC_GA_options2.PopulationSize=50;
PC_GA_options2.MaxGenerations=floor(PC_GA_options.MaxGenerations/10);

if PC_GA_options2.MaxGenerations < 3
    PC_GA_options2.MaxGenerations=3;
end

%1st pass settings
PCFit_start.search=0; %set to zero to trigger this to be a search with an initial seed idea
PCFit_start.PCloc=PCin.start;

%build the outputs
peak_centres=cell(PC_num,1);

%pass 1
PC_out_all=zeros(3,num_phases,PC_num);
PC_out_err=zeros(1,num_phases,PC_num);

%pass 2
PC_out_all2=zeros(3,num_phases,PC_num);
PC_out_err2=zeros(1,num_phases,PC_num);


parfor n=1:PC_num
    %read pattern
    p=PC_patnum(n);
    Peak_Centre_full=Peak_Centres(:,:,p);
    Peak_Centre_ok=Peak_Centre_full(1:Peak_Bands(p),:);
    
    %save into array for the 2nd pass  
    peak_centres{n}=Peak_Centre_ok;
   
    %find the pattern centre
    [PC_out,PC_err]=PC_subfind(PCFit_start,PC_paty,PC_patx,Peak_Centre_ok,outputs,Crystal_LUT,Crystal_UCell,num_phases,n,PC_GA_options,PCin.range);
    
    %store the output
    PC_out_all(:,:,n)=PC_out;
    PC_out_err(1,:,n)=PC_err;
   
end

%pack the data up after the search - important for parfor transparency
PC_data.PC_out_all=PC_out_all;
PC_data.PC_out_err=PC_out_err;
PC_data.PC_patnum=PC_patnum;
PC_data.PC_gridx=PC_gridx;
PC_data.PC_gridy=PC_gridy;
PC_data.PC_patx=PC_patx;
PC_data.PC_paty=PC_paty;


%fit the solved PC data to a grid
[PCFit_out]=pc_fit(PC_num,PC_data,MapArea);

% Run a second pass to refine
%involves an update the pc hunting based upon the first grid fit
parfor n=1:PC_num
    [PC_out,PC_err]=PC_subfind(PCFit_out,PC_paty,PC_patx,peak_centres{n},outputs,Crystal_LUT,Crystal_UCell,num_phases,n,PC_GA_options,PC_subsize);
    PC_out_all2(:,:,n)=PC_out;
    PC_out_err2(1,:,n)=PC_err;   
end

%pack the data back up
PC_data2=PC_data;
PC_data2.PC_out_all=PC_out_all2;
PC_data2.PC_out_err=PC_out_err2;

%perform the second fitting
[PCFit_out2]=pc_fit(PC_num,PC_data2,MapArea);

PCOut.Data_1st=PC_data;
PCOut.Fit_1st=PCFit_out;
PCOut.PC.Data_2nd=PC_data2;
PCOut.Fit_2nd=PCFit_out2;

end

function [PCFit_out]=pc_fit(PC_num,PC_data,MapArea)
%fit the plane
PC_ok=zeros(PC_num,3);
PC_err=zeros(PC_num,1);
PC_phase=zeros(PC_num,1);

for n=1:PC_num
    [mv,ix]=nanmin(PC_data.PC_out_err(:,:,n));
    if ~isnan(mv) && mv < 1E-2 %1E-2 in radian MAD
        PC_ok(n,:)=PC_data.PC_out_all(:,ix,n);
        PC_err(n)=mv;
        PC_phase(n)=ix;
    else
        PC_ok(n,:)=[NaN,NaN,NaN];
        PC_err(n)=NaN;
        PC_phase(n)=0;
    end
end

%remove the dodgy points
PC_fit=PC_ok(PC_phase ~= 0,:);
PC_x=PC_data.PC_patx(PC_phase ~= 0);
PC_y=PC_data.PC_paty(PC_phase ~= 0);

%fit the data to a plane
x_fun=robustfit([PC_x,PC_y],PC_fit(:,1));
y_fun=robustfit([PC_x,PC_y],PC_fit(:,2));
z_fun=robustfit([PC_x,PC_y],PC_fit(:,3));

%calculate the full map data
PCFit_out.PCx_map=MapArea.XBeam_Map*x_fun(2)+x_fun(1);
PCFit_out.PCy_map=MapArea.YBeam_Map*y_fun(3)+y_fun(1);
PCFit_out.PCz_map=MapArea.YBeam_Map*z_fun(3)+z_fun(1);

PCFit_out.search=1;
PCFit_out.precentpoints=length(PC_x(:))/PC_num;


end

function [PC_out,PC_err]=PC_subfind(PCFit_out,PC_paty,PC_patx,Peak_Centre_ok,outputs,Crystal_LUT,Crystal_UCell,num_phases,n,PC_GA_options,sub_refine)
%find the pattern centre using a genetic algorithm 
%for all points in the submap

PC_out=zeros(3,num_phases);
PC_err=zeros(1,num_phases);
    
if PCFit_out.search == 1
    PCx=PCFit_out.PCx_map(PC_paty(n),PC_patx(n));
    PCy=PCFit_out.PCy_map(PC_paty(n),PC_patx(n));
    PCz=PCFit_out.PCz_map(PC_paty(n),PC_patx(n));
    PC_loc=[PCx,PCy,PCz];
else
    PC_loc=[PCFit_out.PCloc];
end

    PC_GA_ub=PC_loc+sub_refine;
    PC_GA_lb=PC_loc-sub_refine;
    
    %find the pattern centre
    for num_P=1:num_phases
        FitFunc = @(PC_test) PC_GAOpt( PC_test,Peak_Centre_ok,outputs.size,Crystal_LUT,Crystal_UCell,num_P);
        [PC_out(:,num_P), PC_err(num_P), ] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
    end
    
end

function [PC_patnum,PC_patx,PC_paty,PC_num,PC_gridx,PC_gridy]=PC_GridBuild(MapArea,PC_array)

%build the grid of locations to search for pattern centre
PC_xline=round(linspace(MapArea.X_bline(1),MapArea.X_bline(end),PC_array(1)));
PC_yline=round(linspace(MapArea.Y_bline(1),MapArea.Y_bline(end),PC_array(2)));
[PC_gridx,PC_gridy]=meshgrid(PC_xline,PC_yline);
PC_num=size(PC_gridx(:),1);
PC_patnum=zeros(PC_num,1);
PC_patx=PC_patnum;
PC_paty=PC_patnum;

for n=1:PC_num
    PC_patnum(n)=MapArea.PMap(PC_gridy(n),PC_gridx(n));
    PC_patx(n)=MapArea.XBeam_Map(PC_gridy(n),PC_gridx(n));
    PC_paty(n)=MapArea.YBeam_Map(PC_gridy(n),PC_gridx(n));
end



end


