function [ EBSP ] = EBSP_Gnom( PatternInfo,PCin )
%EBSP_GNOM Generate the gnomonic coordinate set for a detector geometry
%   Inputs = Decector - a structure with fields
%   
%   using Bruker coordinate systems
%   PC(3)=0.5; %in units of screen height
%   PC(1)=0.5; %in units of screen width, from left of exp. pattern
%   PC(2)=0.5; %in units of screen height, from top of exp. pattern
%   PatternInfo.size = [width height] % width in pixels
%
%   Outputs = EBSP
%   In Gnomonic coordinate system

if isfield(PatternInfo,'ScreenWidth')
    PatternInfo.size(2)=PatternInfo.ScreenWidth;
    PatternInfo.size(1)=PatternInfo.ScreenHeight;
end


if isfield(PCin,'PCX')
    PC(1)=PCin.PCX;
    PC(2)=PCin.PCY;
    PC(3)=PCin.DD;
else
    PC=PCin;
end

EBSP.size=PatternInfo.size;

EBSP.ScreenAspect=double(PatternInfo.size(2)/PatternInfo.size(1));

EBSP.y_gn_max= PC(2)/PC(3);
EBSP.y_gn_min= -(1.0-PC(2))/PC(3);
EBSP.x_gn_max= ((1.0-PC(1))*EBSP.ScreenAspect)/PC(3);
EBSP.x_gn_min= -((PC(1))*EBSP.ScreenAspect)/PC(3);
EBSP.x_screen=linspace(EBSP.x_gn_min,EBSP.x_gn_max,double(PatternInfo.size(2)));
EBSP.y_screen=linspace(EBSP.y_gn_min,EBSP.y_gn_max,double(PatternInfo.size(1)));
EBSP.PC=[PC(1) PC(2) PC(3)];

%set up the screen for interpolation
[EBSP.ypts_screen,EBSP.xpts_screen]=ndgrid(EBSP.y_screen,EBSP.x_screen);
EBSP.r = [EBSP.xpts_screen(:), EBSP.ypts_screen(:), EBSP.ypts_screen(:)*0+1].*1./sqrt((EBSP.xpts_screen(:).^2+EBSP.ypts_screen(:).^2+1));

end

