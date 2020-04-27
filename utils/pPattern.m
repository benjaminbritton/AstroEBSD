function [I1] = pPattern(Pattern,EBSD_Geometry)
%PPATTERN Summary of this function goes here
%   Detailed explanation goes here
if min(size(Pattern)) == 1
    Pattern=reshape(Pattern,EBSD_Geometry.size(1),EBSD_Geometry.size(2));
end

I1=imagesc(EBSD_Geometry.x_screen,EBSD_Geometry.y_screen,Pattern); axis image; axis xy; colormap('gray'); 
end

