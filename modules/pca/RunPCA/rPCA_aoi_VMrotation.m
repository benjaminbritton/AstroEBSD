% rPCA_tpm_VMrotation  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%This generates a structure VMOutput, which contains all the rotated
%patterns, spectra, factor maps.

%For a given area, number of comps to retain, weighting.


%% Rotate the factors to improve contrast in the factors (i.e. make the indexable)
clear VMOutput
VMOutput=struct();

%%
% Fix normalisation issue
Output.coeff1=Output.coeff;
Output.coeff=Output.coeff(:,1:ret_comps);
h = sqrt(sum(Output.coeff.^2, 2));
zerolocs=find(h==0);
Output.coeff(zerolocs,:)=1./double(ret_comps);


if PCA_Setup.cropping==1;
    VMOutput.width=MapInfo.cropped_width;
    VMOutput.height=MapInfo.cropped_height;
end

if PCA_Setup.cropping==0;
    VMOutput.width=MapInfo.MapSize2;
    VMOutput.height=MapInfo.MapSize1;
end

%disp(ret_comps);
VMOutput.NumPCA=ret_comps;

%rotate the factors
[VMOutput.coeffVM, VMOutput.RotVM] = rotatefactors(Output.coeff(:,1:VMOutput.NumPCA),'Normalize','on','Method','varimax', 'Maxit', 10000, 'Reltol', 0.5*sqrt(eps));
VMOutput.scoreVM=Output.score(:,1:VMOutput.NumPCA)*VMOutput.RotVM;

%coeffVM returns coefficients for every component (numPCA) for every map
%pixel. Therefore for numPCA=36 has dimensions (w*h) x 36.
%RotVM is the rotation matrix applied
%scoreVM are our principal components rotated, so the matrix of PC kikuchi
%patterns

%Sort the factor loadings to make them positive

% preallocate the factor maps
VMOutput.PCA_VM_Map_n=zeros([VMOutput.height VMOutput.width VMOutput.NumPCA]);
VMOutput.PCA_VM_num=zeros([VMOutput.height VMOutput.width]);
VMOutput.PCA_VM_val=zeros([VMOutput.height VMOutput.width]);

if PCA_Setup.PCA_EBSD==1;
    VMOutput.PCA_VM_Pat_n=zeros(EBSD_Info.PatSizeW,EBSD_Info.PatSizeH,VMOutput.NumPCA);
end

%sort the loadings
if PCA_Setup.PCA_EBSD==1;
    VMOutput.kikuchi_VM_n=VMOutput.scoreVM(1:(EBSD_Info.PatSizeH.*EBSD_Info.PatSizeW),:);
    if PCA_Setup.PCA_EDX==1;
        VMOutput.spectrum_VM_n=VMOutput.scoreVM((EBSD_Info.PatSizeH.*EBSD_Info.PatSizeW)+1:end,:);
    end
end

if PCA_Setup.PCA_EDX==1;
    if PCA_Setup.PCA_EBSD==0;
        VMOutput.spectrum_VM_n=VMOutput.scoreVM;
    end
end

%% go through and add zeros where patterns weren't loaded in the aoi
VMOutput.coeffVM2=zeros(VMOutput.width*VMOutput.height,VMOutput.NumPCA);
[x_grid,y_grid]=meshgrid(tile.colstart(tiling_no):tile.colfin(tiling_no),tile.rowstart(tiling_no):tile.rowfin(tiling_no));
xn=size(x_grid,2);
yn=size(x_grid,1);

for i=1:size(VMOutput.coeffVM,1)
    
        x=tile.xloaded(i);
        y=tile.yloaded(i);
        ind=sub2ind([VMOutput.height,VMOutput.width],y,x);
    
        VMOutput.coeffVM2(ind,:)=VMOutput.coeffVM(i,:);
        tile.xy_loc(i)=ind;
        
end

%%
%reshaping
for n=1:VMOutput.NumPCA
    
    if PCA_Setup.PCA_EBSD==1; %reshape PCs into PC kikuchi patterns
    VMOutput.PCA_VM_Pat_n(:,:,n)=reshape(VMOutput.kikuchi_VM_n(:,n),[EBSD_Info.PatSizeW,EBSD_Info.PatSizeH]);
    end
    
    % Below is just function of space, so no if statements required
    %how much each pixel in map is loaded by pattern n
    VMOutput.Map_n=reshape(VMOutput.coeffVM2(:,n),[VMOutput.height VMOutput.width]); % map_n is spatial, independent of kikuchis or spectra
    VMOutput.PCA_VM_Map_n(:,:,n)=VMOutput.Map_n;
    
    %PCA_VM_num = map of dominant factor:
    %as we go throuigh n (the PCs) if loading is greater for this n, then
    %take this value. 
    %number = dominant PC; val = actual value of loading of that PC.
    VMOutput.PCA_VM_num(abs(VMOutput.Map_n) > abs(VMOutput.PCA_VM_val)) = n;
    VMOutput.PCA_VM_val(abs(VMOutput.Map_n) > abs(VMOutput.PCA_VM_val)) = VMOutput.Map_n(abs(VMOutput.Map_n) > abs(VMOutput.PCA_VM_val));
    
    VMOutput=rmfield(VMOutput,'Map_n'); %only a function of one PC so tidy away...
end
clear n

%flip the loadings if needed
for n=1:VMOutput.NumPCA
    
    if PCA_Setup.PCA_EBSD==1
        VMOutput.meanV=mean(VMOutput.PCA_VM_val(VMOutput.PCA_VM_num==n));
        if VMOutput.meanV < 0
            VMOutput.PCA_VM_Pat_n(:,:,n)=-VMOutput.PCA_VM_Pat_n(:,:,n);
            VMOutput.PCA_VM_Map_n(:,:,n)=-VMOutput.PCA_VM_Map_n(:,:,n);
            VMOutput.PCA_VM_val(VMOutput.PCA_VM_num==n)=-VMOutput.PCA_VM_val(VMOutput.PCA_VM_num==n);
        end
    end
    
    if PCA_Setup.PCA_EDX==1
        VMOutput.meanSpec=mean(VMOutput.spectrum_VM_n(:,n));
        if VMOutput.meanSpec < 0
            VMOutput.spectrum_VM_n(:,n)=-VMOutput.spectrum_VM_n(:,n);
        end
    end
end

clear n

% PCA radon
%coeff = maps
%scores = patterns

% if PCA_radon == 1
%     pTime('PCA Radon',t1);
%     testRadon2=reshape(testRadon,size(testRadon,1),size(testRadon,2)*size(testRadon,3));
%     [coeffr,scorer,latent,~,explained]=pca(testRadon2, 'Centered',false, 'NumComponents',NumPCA);
%     %rotate the factors
%     [coeffVMr, RotVMr] = rotatefactors(coeffr(:,1:NumPCA)+mean(mean(coeffr)),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
%     coeffVMr=coeffVMr-mean(mean(coeffr));
%     scoreVMr=scorer(:,1:NumPCA)*RotVMr;
% end

% count the number of pixels changed by a median filter
VMOutput.filtered_map=medfilt2(VMOutput.PCA_VM_num);
VMOutput.filterdifference=imbinarize(VMOutput.filtered_map-VMOutput.PCA_VM_num);
VMOutput.filterdifference=size(find(VMOutput.filterdifference),1);

% calculate the ratio of 1st to 2nd strongest component at each pixel,
% averaged over all pixels.
[vals,locs]=max(VMOutput.coeffVM,[],2);
VMOutput.strongestcoeff=mean(vals);

%create dummy coeffVM with maxima removed
VMOutput.coeffVM_1=VMOutput.coeffVM;
VMOutput.coeffVM_1(VMOutput.coeffVM_1==vals)=0;
[vals2,locs2]=max(VMOutput.coeffVM_1,[],2);
% VMOutput.coeffVM_1(bsxfun(@eq, VMOutput.coeffVM_1, max(VMOutput.coeffVM_1,[],2))) = -Inf;
% %secondmax
VMOutput.secondcoeff=mean(vals2);

VMOutput=rmfield(VMOutput,'coeffVM_1');
clear vals locs vals2 locs2 h zerolocs

Output.coeff=Output.coeff1;
