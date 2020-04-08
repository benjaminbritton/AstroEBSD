function [cmap] = violetFire(Nc)

if nargin == 0
    Nc = 256;
end

% Make custom colormap
% [R G B fraction]
c = [0 0 0 0;
    .2 0 .4 .25;
    1 0 0 .5;
    1 .6 0 .65;
    1 1 0 .75;
    1 1 1 1];
cmap = zeros(Nc,3);
for a0 = 1:(size(c,1)-1)
    f1 = round(c(a0,4)*Nc+1);
    f2 = round(c(a0+1,4)*Nc);
    inds = f1:f2;
    cnew = [linspace(c(a0,1),c(a0+1,1),length(inds))' ...
        linspace(c(a0,2),c(a0+1,2),length(inds))' ...
        linspace(c(a0,3),c(a0+1,3),length(inds))'];
    cmap(inds,:) = cnew;
end
% smooth colormap in hsv space?
% cmap = rgb2hsv(cmap);
% cmap(:,1) = smooth(cmap(:,1)-.5,11,'moving');
% cmap = hsv2rgb(cmap);
cmap(:,1) = smooth(cmap(:,1),11,'moving');
cmap(:,2) = smooth(cmap(:,2),11,'moving');
cmap(:,3) = smooth(cmap(:,3),11,'moving');
cmap(cmap<0) = 0;
cmap(cmap>1) = 1;

end