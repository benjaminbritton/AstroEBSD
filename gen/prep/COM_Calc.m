function [output] = COM_Calc(image1)
%COM_CALC Calculate the centre of mass of an image
mass = sum(image1(:));
[yy,xx] = ndgrid(1:size(image1,1),1:size(image1,2));
y = sum(yy(:).*image1(:))/mass;
x = sum(xx(:).*image1(:))/mass;
output = [mass,x,y];
end

