function [angle_matrix,xgrid,ygrid] = TKD_scattergrid(PCx,PCy,DD,npix)
%TKD Scattering Grid
%calculates the scattered grid pixel positions

[xgrid,ygrid] = meshgrid(1:npix,1:npix);
angle_matrix = atan(sqrt((xgrid - PCx).^2 + (ygrid - PCy).^2) * 55 / 1000 ./ DD) / pi * 180;

end

