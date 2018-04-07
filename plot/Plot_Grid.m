function PlotDetails=Plot_Grid(PlotSet,MapOut,ScaleDetails)
%plot data in a tidy grid image

% num_grid=PlotSet.Grid(1)*PlotSet.Grid(2); %Y, X

psize=[PlotSet.Grid(1)*3/5*10,PlotSet.Grid(2)*10];
f1=figure('PaperSize',psize*3/4,'Name',PlotSet.Title,'Clipping','off','Visible','off');
set(gcf,'Units','centimeters','Position',[5   3   psize(2)    psize(1)]);

% PlotSet.Grid=[3,2];
% PlotSet.Data={MapOut.IQ,MapOut.BQ,MapOut.BandNumber,MapOut.Err*180/pi,MapOut.Phase};
% PlotSet.Locations=[1,2,3,4,5];
% PlotSet.Title='Quality Data';
% PlotSet.FigTitles={'Image Quality','Band Slope','Band Number','Mean Angular Error','Phase'};
% PlotSet.CRanges={1,1,1,[0 3],0};
% PlotSet.ImageType=[1,1,1,1,1]; %1 = imagesc, 2 = image
% PlotSet.CMap='grey';



for n=1:size(PlotSet.Data,2)
    s1(n)=subplot(PlotSet.Grid(1),PlotSet.Grid(2),PlotSet.Locations(n));
    switch PlotSet.ImageType(n)
        case 1 %imagesc
            i1(n)=imagesc(MapOut.X_axis,MapOut.Y_axis,PlotSet.Data{n});
            if PlotSet.CRanges{n} == 0
                c1=colorbar;
                c1.Visible='off';
            elseif PlotSet.CRanges{n} == 1
                c1=colorbar;
            else
                c1=colorbar;
                caxis(PlotSet.CRanges{n});
            end
        case 2 %image
            i1(n)=image(MapOut.X_axis,MapOut.Y_axis,PlotSet.Data{n});
            c1=colorbar;
            c1.Visible='off';
    end
    
    if isfield(PlotSet,'Alphas')
        if ~isempty(PlotSet.Alphas{n})
            i1(n).AlphaData=PlotSet.Alphas{n};
        end
    end
    
    axis equal; axis tight; axis off;
    title(PlotSet.FigTitles{n});
    hold on; 
    plot(ScaleDetails.X,ScaleDetails.Y,'w','LineWidth',3);
    plot(ScaleDetails.X,ScaleDetails.Y,'k','LineWidth',1);
    
    if strcmpi(MapOut.CoordSystems,'TRZP')
        s1(n).XDir='reverse';
    end
end
if isfield(PlotSet,'CMap')
colormap(PlotSet.CMap);
end

a1=annotation(f1,'textbox',[0.1 0.04 0.2 0.04035],'String',ScaleDetails.Text,...
    'FitBoxToText','off','FontSize',10,'EdgeColor','none','Position',[0.1 0.04 0.5 0.04],'FitBoxToText','on');
switch PlotSet.Grid(1)
    case 1
    a1.Position=[0.1 0.05 0.5 0.04];
    case 3
    a1.Position=[0.1 0.05 0.5 0.04];
end

f1.Visible='on';

PlotDetails.images=i1;
PlotDetails.subplots=s1;
PlotDetails.annotations=a1;
PlotDetails.figures=f1;
