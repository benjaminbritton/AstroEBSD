function [PatternData] = fSelectPoints(MapInfo,PhaseInput,no_points)
%FSELECTPOINTS Summary of this function goes here
%   Detailed explanation goes here

disp('Please select points from the GUI window');
Radon=MapInfo.Data_InputMap.RadonQuality;

figure
I=imagesc(Radon);
colormap('gray')

%t=text(8,10,['Select ',num2str(no_points),' points for refinement...'],'Color','blue','FontSize',14);
t=annotation('textbox',[0.15,0.4,1,0.5],'String',['Select ',num2str(no_points),' ',PhaseInput,' points for refinement...'],'FitBoxToText','on');
t.Color=[1,1,1];
t.BackgroundColor=[0,0,0];
t.FontName='Arial';
t.FontSize=13;

t=annotation('textbox',[0.15,0.1,1,0.1],'String',['For best performance, evenly spaced and > 1 orientation.'],'FitBoxToText','on');
t.Color=[1,1,1];
t.BackgroundColor=[0,0,0];
t.FontName='Arial';
t.FontSize=9;

drawnow

for i = 1:no_points
    p(i) = images.roi.Point;
end

for i = 1:no_points
    draw(p(i))
end

%convert to useful data to output
for i =1:no_points
    locations(i,:)=floor(p(1,i).Position);
    
    %Bruker Data
    PatternData.PC_start(i,:)=[MapInfo.Data_InputMap.PCX(locations(i,2),locations(i,1)),MapInfo.Data_InputMap.PCY(locations(i,2),locations(i,1)),MapInfo.Data_InputMap.DD(locations(i,2),locations(i,1))]; %initial value for PC
    PatternData.Eulers(i,:)=[MapInfo.Data_InputMap.phi1(locations(i,2),locations(i,1))*degree,MapInfo.Data_InputMap.PHI(locations(i,2),locations(i,1))*degree,MapInfo.Data_InputMap.phi2(locations(i,2),locations(i,1))*degree]; %initial value for Eulers
    
    %Pattern Number
    PatternData.P(i)=MapInfo.Data_InputMap.PMap(locations(i,2),locations(i,1));
    
end

%disp('Point selection complete.');

end

