%% Plot some data
%error filter

if Settings_Mode == 4 %full map
    Plots_FilterMap=Plot_FilterGen(Data_OutputMap,Settings_PlotFilters);
    %build the scale bar
    [Plots_Scale]=Plot_SetScaleBar(Data_OutputMap);
    
    % Quality
    PlotSetQuality.Title='Quality Data';
    PlotSetQuality.Grid=[3,2];
    PlotSetQuality.Data={Data_OutputMap.IQ,Data_OutputMap.BQ,Data_OutputMap.BandNumber,Data_OutputMap.Err*180/pi,Data_OutputMap.Phase};
    PlotSetQuality.Locations=[1,2,3,4,5];
    PlotSetQuality.FigTitles={'Image Quality','Band Slope','Band Number','Mean Angular Error','Phase'};
    PlotSetQuality.CRanges={1,1,1,[0 3],0};
    PlotSetQuality.ImageType=[1,1,1,1,1]; %1 = imagesc, 2 = image
    PlotSetQuality.CMap='gray';
    Plots_Quality=Plot_Grid(PlotSetQuality,Data_OutputMap,Plots_Scale);
    
    print(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_Quality']),'-dpng','-r300')
    savefig(Plots_Quality.figures,fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_Quality']));
    
    % IPF
    PlotSetIPF.Title='Inverse Pole Figures Data';
    PlotSetIPF.Grid=[1,3];
    PlotSetIPF.Data={Plot_IPF.X,Plot_IPF.Y,Plot_IPF.Z};
    PlotSetIPF.Locations=[1,2,3];
    PlotSetIPF.FigTitles={'IPF-X','IPF-Y','IPF-Z'};
    PlotSetIPF.ImageType=[2,2,2]; %1 = imagesc, 2 = image
    PlotSetIPF.Alphas={Plots_FilterMap,Plots_FilterMap,Plots_FilterMap};
    Plots_IPF=Plot_Grid(PlotSetIPF,Data_OutputMap,Plots_Scale);
    
    print(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_IPF']),'-dpng','-r300')
    savefig(Plots_IPF.figures,fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_IPF']));
    
    % PC
    PlotSetPC.Title='Beam Position & Pattern Centre';
    PlotSetPC.Grid=[3,3];
    PlotSetPC.Data={Data_OutputMap.X_co,Data_OutputMap.Y_co,Data_OutputMap.P_co,PCOut.Fit_2nd.PCx_map,...
        PCOut.Fit_2nd.PCy_map,PCOut.Fit_2nd.PCz_map,...
        Data_OutputMap.XGrid,Data_OutputMap.YGrid};
    PlotSetPC.Locations=[1,2,3,4,5,6,7,8];
    PlotSetPC.CRanges={1,1,1,1,1,1,1,1};
    PlotSetPC.FigTitles={'Beam X','Beam Y','Pattern Number','PC_x','PC_y','DD','Map X','Map Y'};
    PlotSetPC.ImageType=[1,1,1,1,1,1,1,1];
    Plots_PC=Plot_Grid(PlotSetPC,Data_OutputMap,Plots_Scale);
    
    print(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_PC']),'-dpng','-r300')
    savefig(Plots_Quality.figures,fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_PC']));
end

%% One pattern plot
EBSP_OneFigure=Plot_SinglePattern(EBSP_One,Crystal_UCell,Crystal_LUT,InputUser.PatternPhase);
print(fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_Pattern']),'-dpng','-r300')
savefig(EBSP_OneFigure.figure,fullfile(InputUser.FolderOut,[InputUser.EBSD_File '_' InputUser.FileOut '_Pattern']));



