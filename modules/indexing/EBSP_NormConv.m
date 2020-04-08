function [ nhat_gnom] = EBSP_NormConv( Peak_Centre_ok,ScreenSizeS,PC)
%EBSP_NORMCONV Summary of this function goes here
%   Detailed explanation goes here

%set up the screen coordinate systems in gnomonic
ScreenAspect=ScreenSizeS(2)/ScreenSizeS(1);
Gnom_Step=(1/PC(3))/(ScreenSizeS(1)-1);

IC_Gnomy=((floor((ScreenSizeS(1)+1)/2)-1)/(ScreenSizeS(1)-1)-(1.0-PC(2)))/PC(3);
IC_Gnomx=ScreenAspect*((floor((ScreenSizeS(2)+1)/2)-1)/(ScreenSizeS(2)-1)-PC(1))/PC(3);

c_theta=cos(Peak_Centre_ok(:,1)*pi/180);
t_theta=tan(Peak_Centre_ok(:,1)*pi/180);

va1=[zeros(size(Peak_Centre_ok,1),1),-IC_Gnomx+IC_Gnomy.*t_theta-Gnom_Step*Peak_Centre_ok(:,2)./c_theta,t_theta];
va3=[t_theta,ones(size(Peak_Centre_ok,1),1),zeros(size(Peak_Centre_ok,1),1)];



% va3=va3./repmat(sqrt(dot(va3,va3,2)),1,3);
% va1=va1./repmat(sqrt(dot(va1,va1,2)),1,3);

%calculate the unit normal vectors
nhat_gnom=cross(va1,va3,2);
nhat_gnom=nhat_gnom./repmat(sqrt(dot(nhat_gnom,nhat_gnom,2)),1,3);

% nhat_2=[ones(size(Peak_Centre_ok,1),1),t_theta,-IC_Gnomx+IC_Gnomy.*t_theta-Gnom_Step*Peak_Centre_ok(:,2)./c_theta];
% nhat_2=nhat_2./repmat(sqrt(dot(nhat_2,nhat_2,2)),1,3);

end

