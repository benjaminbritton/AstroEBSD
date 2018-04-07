function EBSP_OneFigure=Plot_SinglePattern(EBSP_One,Crystal_UCell,Crystal_LUT,Num_Phase)

PlotSet.Grid=[2 2];
psize=[PlotSet.Grid(1)*3/5*10,PlotSet.Grid(2)*10];
EBSP_OneFigure.figure=figure('PaperSize',psize*3/4,'Name','Single Pattern Data','Clipping','off','Visible','off');
set(gcf,'Units','centimeters','Position',[5   3   30    18]);

EBSP_OneFigure.s1(1)=subplot(2,2,1);
imagesc(EBSP_One.PatternIn); axis equal; colormap('gray'); axis xy; axis equal; axis tight;
title('Raw Input');

EBSP_OneFigure.s1(2)=subplot(2,2,2);
imagesc(EBSP_One.R_theta,EBSP_One.R_rho,EBSP_One.R_EBSP); colormap('gray'); axis xy; axis equal; axis tight;
hold on; scatter(EBSP_One.Peak_Centre(:,1),EBSP_One.Peak_Centre(:,2),'r');
xlabel('Theta ( \circ )'); ylabel('Rho');
title('Radon Transform with peaks found');

EBSP_OneFigure.s1(3)=subplot(2,2,3);
Plot_EBSPNorms( EBSP_One.PatternCor,EBSP_One.PatternGeometry,EBSP_One.nhat_gnom,EBSP_OneFigure.s1(3));
title('Background corrected & RT Bands'); 
axis image; 
xlim([EBSP_One.PatternGeometry.x_gn_min EBSP_One.PatternGeometry.x_gn_max])
ylim([EBSP_One.PatternGeometry.y_gn_min EBSP_One.PatternGeometry.y_gn_max])
axis xy;

EBSP_OneFigure.s1(4)=subplot(2,2,4);
Plot_EBSPAnnotated( EBSP_One.PatternCor,EBSP_One.PatternGeometry,EBSP_One.nhat_gnom,EBSP_One.rotdata{Num_Phase}.detector,Crystal_UCell{Num_Phase},Crystal_LUT{Num_Phase}.family_norm_list,EBSP_OneFigure.s1(4));
title('Indexed Pattern');
axis image; 
xlim([EBSP_One.PatternGeometry.x_gn_min EBSP_One.PatternGeometry.x_gn_max])
ylim([EBSP_One.PatternGeometry.y_gn_min EBSP_One.PatternGeometry.y_gn_max])
axis xy;

EBSP_OneFigure.figure.Visible='on';