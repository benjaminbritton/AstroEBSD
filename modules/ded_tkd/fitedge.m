function [fitlinex, fitliney] = fitedge(x1,x2,y1,y2)
fitslope = (y2 - y1) / (x2 - x1);

fitlinex = linspace(x1,x2,(abs(x2 - x1) + 1));
fitliney = y1 + fitslope .* (fitlinex - x1);
end
