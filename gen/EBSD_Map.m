function [AreaData] = EBSD_Map( MapData,MicroscopeData )
%EBSD_Map Generate pattern maps that link pattern number to position
%INPUTS
%MapData        = MapData structure
%MicrospeData   = Microscope information
%
%Outputs
%AreaData = structure with data in map form
%
%Assumes map is rectangular

NCOLS=double(MicroscopeData.NCOLS);
NROWS=double(MicroscopeData.NROWS);
max_pats=NCOLS*NROWS;

PatList2Map=@(data,col,row)fliplr(rot90(reshape(double(data),[double(col) double(row)]),3));

plist=1:max_pats;
plist=plist(:);

XBeam_Map=PatList2Map(MapData.XBeam,NCOLS,NROWS);
YBeam_Map=PatList2Map(MapData.YBeam,NCOLS,NROWS);
PMap=PatList2Map(plist,NCOLS,NROWS);

X_axis=XBeam_Map(1,:); X_axis=X_axis(:);
Y_axis=YBeam_Map(:,1); Y_axis=Y_axis(:);

AreaData.XStep=double(round(MicroscopeData.XSTEP/1E-4)*1E-4);
AreaData.YStep=double(round(MicroscopeData.YSTEP/1E-4)*1E-4);

AreaData.X_bline=X_axis;
AreaData.Y_bline=Y_axis;
AreaData.X_axis=X_axis(:).*AreaData.XStep;

AreaData.CoordSystems=MicroscopeData.CoordSystems;


if strcmpi(MicroscopeData.CoordSystems,'TRZP') == 1 %top right z plus
    AreaData.X_axis=flipud(AreaData.X_axis);
end
AreaData.Y_axis=Y_axis(:).*AreaData.YStep;

AreaData.XBeam_Map=XBeam_Map;
AreaData.YBeam_Map=YBeam_Map;
AreaData.PMap=PMap;
AreaData.max_pats=max_pats;
AreaData.xpts=size(X_axis,1);
AreaData.ypts=size(Y_axis,1);

end

