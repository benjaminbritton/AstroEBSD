cd(printingfol)

figure;
imagesc(Plots.SEMImageCropped(:,:,tiling_no)); axis image off; colormap('gray'); title('SEM Image');
print(gcf,[name,'SEM_cropped'],'-dpng','-r300');

figure
imagesc(Plots.RadonCropped(:,:,tiling_no));axis image off;colormap('gray'); title('IQ');
print(gcf,[name,'_PQcropped'],'-dpng','-r300');

% if PCA_Setup.PCA_EDX==1;
%     figure
%     imagesc(Plots.EDSInteg_cropped)
%     print(gcf,[name,'EDS_integ_cropped'],'-dpng','-r300'); 
% end
    

% %% Plot rot iterations
% if exist('RotIterations_par','var')==1
% 
%     figure
%     yyaxis left
%     p1=plot(RotIterations_par.NumPCA_it, RotIterations_par.ratio);
%     hold on
%     xlabel('Number of components retained');
%     xlim([min(RotIterations_par.NumPCA_it),max(RotIterations_par.NumPCA_it)]);
%     ylabel('Ratio of 1st to 2nd component')
%     yyaxis right
%     p2=plot(RotIterations_par.NumPCA_it,RotIterations_par.strongest_it);
%     hold on
%     p3=plot(RotIterations_par.NumPCA_it, RotIterations_par.secondstrongest_it);
%     ylabel('Coefficient')
%     legend([p1, p2, p3],'Ratio','Strongest', 'Second strongest','Location','SouthEast')
% 
%     print(gcf,[name,'CompRatios_',num2str(RotIterations_par.ret_comp_range(1)),'-',num2str(RotIterations_par.ret_comp_range(2))],'-dpng','-r300');
%     
%     figure
%     plot(RotIterations_par.NumPCA_it, RotIterations_par.filtercount);
%     xlabel('Number of components retained')
%     ylabel('Number of pixels removed by median filter')
%     
%     print(gcf,[name,'MedFilterRemoved_',num2str(RotIterations_par.ret_comp_range(1)),'-',num2str(RotIterations_par.ret_comp_range(2))],'-dpng','-r300');
% end
%     
% clear p1 p2 p3
% close all