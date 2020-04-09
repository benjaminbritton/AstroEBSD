function [rotdata,banddata] = Map_Index( AreaData,Peak_Centres,Peak_Bands,num_Phases,PCFit_out,Crystal_UCell,Crystal_LUT,Settings_LUT,PaternInfo,MicroscopeData)
%MAP_INDEX Index a series of theta,rho information for a number of phases
%
%INPUTS
%AreaData = data on the EBSD map and how pattern numbers work
%Peak_Centres = peak centre information - from radon transform
%num_Phases = number of phases to index
%Crystal_UCell = Unit cell information
%Crystal_LUT = Look up table for indexing
%Settings_LUT = Look up table settings
%
%OUTPUTS
%rotdata = data about the indexed case (as cell)
%banddata = band information for indexed case (as cell)

%% Versioning
% v1 = TBB 14/04/2017 - documented code + restructured

%%

%allocate the output cells
banddata=cell(AreaData.max_pats,num_Phases);
rotdata=cell(AreaData.max_pats,num_Phases);

%set up the detector tilt info
R_x=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
rot_det=R_x(MicroscopeData.TotalTilt);

%find the location of the pattern & sort out the pattern centre data
%set up to be suitable for transparency + slicing requirements

y_co=zeros(AreaData.max_pats,1);
x_co=zeros(AreaData.max_pats,1);
PC_pat=zeros(AreaData.max_pats,3);

for n=1:AreaData.max_pats
    [y_co(n),x_co(n)]=find(AreaData.PMap == n);
    PC_pat(n,:)=[PCFit_out.PCx_map(y_co(n),x_co(n)),PCFit_out.PCy_map(y_co(n),x_co(n)),PCFit_out.PCz_map(y_co(n),x_co(n))];
end

parfor n=1:AreaData.max_pats
    Peak_Centre_full=Peak_Centres(:,:,n);
    Peak_Centre_ok=Peak_Centre_full(1:Peak_Bands(n),:);
    
    [ nhat_gnom] = EBSP_NormConv( Peak_Centre_ok,[PaternInfo.size],PC_pat(n,:));
    
    %index for all phases
    for num_P=1:num_Phases
        [rotdata{n,num_P},banddata{n,num_P}]=EBSP_Index(nhat_gnom,Crystal_LUT{num_P},Settings_LUT{num_P}.thresh_trig,Crystal_UCell{num_P},rot_det); %#ok<PFBNS>
    end
end


end

