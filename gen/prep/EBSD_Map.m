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

%update 24/08/2018 (TBB)
% NCOLS=double(MicroscopeData.NCOLS);
%deal with reduced area data - assume rectangular croppping
NCOLS=double(max(MapData.XBeam)-min(MapData.XBeam))+1;
NROWS=double(max(MapData.YBeam)-min(MapData.YBeam))+1;
% NROWS=double(MicroscopeData.NROWS);
max_pats=NCOLS*NROWS;
% if max_pats ~= double(size(MapData.DD,1))
%     error('This map is not rectangular')
% end

minrow=min(MapData.YBeam);
mincol=min(MapData.XBeam);

%% sort the data

fullarray=zeros(MicroscopeData.NROWS,MicroscopeData.NCOLS);
fullarray=fullarray.*nan;

pv_list=sub2ind(double([MicroscopeData.NCOLS MicroscopeData.NROWS]),double(MapData.XBeam),double(MapData.YBeam));
[pv_list2,ix]=sort(pv_list);
pv_list_set=1:size(pv_list2);

AreaData.PMap=reshape(pv_list_set(ix),[NCOLS NROWS])';
   
contents={'DD','MAD','Phase','MADPhase','NIndexedBands','PCX','PCY','PHI','phi2','XBeam','YBeam','phi1','RadonQuality','XSample','YSample'};
for i=1:length(contents)
    area=fullarray;
    try
        if strcmp(contents{i},'PMap')
            vals=pv_list_set(ix);
        else
            vals=MapData.(contents{i});
        end
    catch
       warning([contents{i},' information not loaded!']);
       continue
    end
        
    area(pv_list)=vals;
    area=reshape(area,size(area,2),size(area,1));
    area=area(mincol:mincol+NCOLS-1,minrow:minrow+NROWS-1)';

    if strcmp(contents{i},'XBeam')
        name={'XBeam_Map'};
    elseif strcmp(contents{i},'YBeam')
        name={'YBeam_Map'};
    else 
        name=contents(i);
    end
    AreaData.(name{:})=area;
            
end

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