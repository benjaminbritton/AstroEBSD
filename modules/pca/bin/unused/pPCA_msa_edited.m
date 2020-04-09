%% Plot the PCAs
if plot_pca == 1
    figure;
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(reshape(score(:,n),[PatSizeW,PatSizeH])); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
    end
    
    figure
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(reshape(coeff(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;title(int2str(n));
    end
    
    %plot these new factors
    figure;
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(PCA_VM_Pat_n(:,:,n)); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
    end
    
    figure
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(PCA_VM_Map_n(:,:,n)); axis image; axis ij; axis tight; axis off;title(int2str(n));
    end
    
    %plot the factor map
    f1=figure;
    imagesc(PCA_VM_num); axis image; colormap('colorcube');
end
% 
% if PCA_figures==1
%     for n=1:NumPCA
%     PCPat=PCA_VM_Pat_n(:,:,n);
%     PCPat=PCPat-min(min(PCPat));
%     PCPat=PCPat./max(PCPat(:));
%     imwrite(PCPat,'test1.tif');
%     end
% end
