


%% plot with MTEX

%template data rotation
Template_Phi1=reshape(TemData(1,:),Data_InputMap.ypts,Data_InputMap.xpts); %in rads
Template_PHI=reshape(TemData(2,:),Data_InputMap.ypts,Data_InputMap.xpts);
Template_Phi2=reshape(TemData(3,:),Data_InputMap.ypts,Data_InputMap.xpts);

rot_Template=rotation('Euler',Template_Phi1,Template_PHI,Template_Phi2,'ZXZ');

Template_Quality=reshape(TemData(7,:),Data_InputMap.ypts,Data_InputMap.xpts);

%set MTEX properties
propsT.x=Data_InputMap.XSample;
propsT.y=Data_InputMap.YSample;
propsT.quality=Template_Quality;
phasesT=propsT.quality*0;

EBSD_template=EBSD(rot_Template,phasesT,cs,'options',propsT);

%% Plot the data
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

oM=ipfTSLKey(EBSD_template(cs.mineral));

figure

sp1=subplot(1,3,1);
oM.inversePoleFigureDirection = xvector;
plot(EBSD_template(cs.mineral), oM.orientation2color(EBSD_template(cs.mineral).orientations),'parent',sp1)
title('IPF - X')

sp2=subplot(1,3,2);
oM.inversePoleFigureDirection = yvector;
plot(EBSD_template(cs.mineral), oM.orientation2color(EBSD_template(cs.mineral).orientations),'parent',sp2)
title('IPF - Y')

sp3=subplot(1,3,3);
oM.inversePoleFigureDirection = zvector;
plot(EBSD_template(cs.mineral), oM.orientation2color(EBSD_template(cs.mineral).orientations),'parent',sp3)
title('IPF - Z')
% 
% plot the peak heigght (quality) map
figure
plot(EBSD_template,propsT.quality); title('FFT - pattern quality')
% caxis([0.78 0.84]);
colorbar

