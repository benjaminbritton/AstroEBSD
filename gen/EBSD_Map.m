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
%
%Updated 18/06/2018 (TBB) to deal with the pattern structure more properly

NCOLS=double(MicroscopeData.NCOLS);
NROWS=double(MicroscopeData.NROWS);
max_pats=NCOLS*NROWS;
%%
pv_list=sub2ind(double([MicroscopeData.NCOLS MicroscopeData.NROWS]),double(MapData.XBeam),double(MapData.YBeam));
[pv_list2,ix]=sort(pv_list);
AreaData.PMap=reshape(pv_list2,[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.DD=reshape(double(MapData.DD(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.MAD=reshape(double(MapData.MAD(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.MADPhase=reshape(double(MapData.MADPhase(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.NIndexedBands=reshape(double(MapData.NIndexedBands(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.PCX=reshape(double(MapData.PCX(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.PCY=reshape(double(MapData.PCY(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.PHI=reshape(double(MapData.PHI(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.phi2=reshape(double(MapData.phi2(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.phi1=reshape(double(MapData.phi1(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.RadonQuality=reshape(double(MapData.RadonQuality(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';

AreaData.XBeam_Map=reshape(double(MapData.XBeam(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.YBeam_Map=reshape(double(MapData.YBeam(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.XSample=reshape(double(MapData.XSample(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';
AreaData.YSample=reshape(double(MapData.YSample(ix)),[MicroscopeData.NCOLS MicroscopeData.NROWS])';

% PatList2Map=@(data,col,row)fliplr(rot90(reshape(double(data),[double(col) double(row)]),3));
% 
% plist=1:max_pats;
% plist=plist(:);
% 
% XBeam_Map=PatList2Map(MapData.XBeam,NCOLS,NROWS);
% YBeam_Map=PatList2Map(MapData.YBeam,NCOLS,NROWS);
% PMap=PatList2Map(plist,NCOLS,NROWS);

X_axis=AreaData.XBeam_Map(1,:); X_axis=X_axis(:);
Y_axis=AreaData.YBeam_Map(:,1); Y_axis=Y_axis(:);

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

% AreaData.XBeam_Map=XBeam_Map;
% AreaData.YBeam_Map=YBeam_Map;
% AreaData.PMap=PMap;
AreaData.max_pats=max_pats;
AreaData.xpts=size(X_axis,1);
AreaData.ypts=size(Y_axis,1);

end

