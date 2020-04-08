function [MapOut] = Map_Generate(AreaData,rotdata,Peak_Quality,Crystal_UCell,num_Phases)
%MAP_QUALITYGEN Summary of this function goes here
%www.expmicromech.com
%
%INPUTS
%AreaData = Map information
%rotat data = structure with indexed information
%Peak_Quality = structure with quality data
%num_Phases = number of phases in the map
%
%OUTPUTS
%Rotation_MapS = [3x3,y,x] rotation matrix in the sample frame
%MapOut = structure with mapped outputs
%   IQ = Hough peak height
%   BQ = Band slope
%   Err = Mean angular Error
%   BandNumber = number of bands used for indexing


%% Versioning
%v1 - TBB 14/04/2017
%v1.1 TBB 14/04/2017 
%         = update to structure output, with all the data contained inside
%%

MapOut.IQ=zeros(AreaData.ypts,AreaData.xpts);
MapOut.BQ=MapOut.IQ;
MapOut.Err=MapOut.IQ;
MapOut.BandNumber=MapOut.IQ;

MapOut.phi1=MapOut.IQ;
MapOut.PHI=MapOut.IQ;
MapOut.phi2=MapOut.IQ;

MapOut.GSample=zeros([3,3,size(AreaData.PMap)]);
MapOut.GDetector=MapOut.GSample;
MapOut.X_co=MapOut.IQ;
MapOut.Y_co=MapOut.IQ;
MapOut.P_co=MapOut.IQ;

for n=1:AreaData.max_pats
    [y_co,x_co]=find(AreaData.PMap == n);
    
    MapOut.X_co(y_co,x_co)=x_co;
    MapOut.Y_co(y_co,x_co)=y_co;
    MapOut.P_co(y_co,x_co)=n;
       
    err_list=zeros(num_Phases,1);
    for num_p=1:num_Phases
        err_list(num_p)=rotdata{n,num_p}.error;
    end
    [err,phase]=nanmax(err_list);
    
    MapOut.Phase(y_co,x_co)=phase;
    MapOut.Err(y_co,x_co)=err;
    
    MapOut.IQ(y_co,x_co)=Peak_Quality(n,1);
    MapOut.BQ(y_co,x_co)=Peak_Quality(n,2);

    MapOut.BandNumber(y_co,x_co)=rotdata{n,1}.maxok;

    if phase > 0
        MapOut.GSample(:,:,y_co,x_co)=rotdata{n,phase}.sample;
        MapOut.GDetector(:,:,y_co,x_co)=rotdata{n,phase}.detector;
        MapOut.Crystal(y_co,x_co)=Crystal_UCell{phase}.crystal;
        MapOut.phi1(y_co,x_co)=rotdata{n,phase}.eang(1);
        MapOut.PHI(y_co,x_co)=rotdata{n,phase}.eang(2);
        MapOut.phi2(y_co,x_co)=rotdata{n,phase}.eang(3);
    end

end

%copy a few outputs to make MapOut more useful
MapOut.xpts=max(MapOut.X_co(:));
MapOut.ypts=max(MapOut.Y_co(:));
MapOut.npts=AreaData.max_pats;
MapOut.X_axis=AreaData.X_axis;
MapOut.Y_axis=AreaData.Y_axis;
MapOut.XStep=AreaData.XStep;
MapOut.CoordSystems=AreaData.CoordSystems;

[MapOut.XGrid,MapOut.YGrid]=meshgrid(MapOut.X_axis,MapOut.Y_axis);

end

