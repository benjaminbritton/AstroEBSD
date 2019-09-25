function [DataMap] = bMapSort(MapData,MicroscopeData,DataLine)
%BMAPSORT Sort the data into 2D Maps

%we are going to create a map list of data, based upon the XBeam and YBeam maps
NCOLS=double(max(MapData.XBeam)-min(MapData.XBeam))+1;
NROWS=double(max(MapData.YBeam)-min(MapData.YBeam))+1;

pv_list=sub2ind(double([MicroscopeData.NCOLS MicroscopeData.NROWS]),double(MapData.XBeam),double(MapData.YBeam));
[~,ix]=sort(pv_list);
DataMap=reshape(double(DataLine(ix)),[NCOLS NROWS])';

end

