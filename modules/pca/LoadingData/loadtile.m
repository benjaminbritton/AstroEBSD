% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [Data,tile]=loadtile(tile,tiling_no,EBSD_Info,MapInfo,PCA_Setup,Settings_Cor,t1)
%load a tile of data

% n=PCA_Setup.crop_factor;

tile.xtile=[1:1:PCA_Setup.crop_factor];
tile.ytile=repmat(tile.xtile,1,PCA_Setup.crop_factor);
tile.xtile=repmat(tile.xtile,PCA_Setup.crop_factor,1);
tile.xtile=reshape(tile.xtile,1,PCA_Setup.crop_factor.^2);

%start and end points in row and column
tile.rowstart(tiling_no)=(tile.ytile(tiling_no)-1).*MapInfo.cropped_height+1;
tile.rowfin(tiling_no)=tile.rowstart(tiling_no)+MapInfo.cropped_height-1;
tile.colstart(tiling_no)=(tile.xtile(tiling_no)-1).*MapInfo.cropped_width+1;
tile.colfin(tiling_no)=tile.colstart(tiling_no)+MapInfo.cropped_width-1;

%output as vector EBSP, yi and xi = row and column of the tile
Data_InputMap=MapInfo.Data_InputMap;

%createa list of the chunks to load
[x_grid,y_grid]=meshgrid(tile.colstart(tiling_no):tile.colfin(tiling_no),tile.rowstart(tiling_no):tile.rowfin(tiling_no));
xn=size(x_grid,2);
yn=size(x_grid,1);

%some of the pattern numbers may be nan - ie where the beam was not scanned
%convert to p_numbers
pmap=Data_InputMap.PMap(tile.rowstart(tiling_no):tile.rowfin(tiling_no),tile.colstart(tiling_no):tile.colfin(tiling_no));
p_nums=pmap(~isnan(pmap));

%grab the non-nan elements of pmap
%p_nums=intersect(pmap,p_nums);

p_tot=numel(p_nums);

if PCA_Setup.PCA_EBSD==1
    %read the 1st pattern to get dimensions
    [ pat_example ] = bReadEBSP(MapInfo.EBSPData,1);
    [ pat_example_bg] = EBSP_BGCor( pat_example,Settings_Cor );
    pwid=size(pat_example,2);
    phigh=size(pat_example,1);
end

if PCA_Setup.SpatialKernel==1
    pat_example_bg=pat_example; %lazy coding from me - this is because we want to skip bg-correction until after spatial kernel if required
end

%create a list of pattern numbers that correspond to
%the grid positions you want
% p_nums=zeros(p_tot,1);
% for p=1:p_tot
%     p_nums(p)=Data_InputMap.PMap(y_grid(p),x_grid(p));
% %     y_p(p)=y_grid(p);
% %     x_p(p)=x_grid(p);
% end

%find the continguous blocks
[p_sort,p_order]=sort(p_nums);
p_dif=diff(p_sort);
block_ends=find(p_dif~=1);

if isempty(block_ends)
    block_ends=numel(p_sort);
else
    block_ends(end+1)=numel(p_sort);
end
block_ends=[0;block_ends];

%% Do the EBSD work
if PCA_Setup.PCA_EBSD==1
    %create the empty tile - full pattern resolution
    ebsp_tile=zeros(phigh,pwid,p_tot);
    ebspBG_tile=zeros(size(pat_example_bg,1),size(pat_example_bg,2));

    %read the blocks
    for n=1:numel(block_ends)-1
        [ ebsp_tile(:,:,[block_ends(n)+1]:block_ends(n+1))] = bReadEBSP(MapInfo.EBSPData,[(p_sort(block_ends(n)+1)) (block_ends(n+1)-block_ends(n))]);
    end

    %background correct
    if PCA_Setup.SpatialKernel==0 %don't bg correct YET if going to use a spatial kernel later
        for n=1:p_tot
            [ ebspBG_tile(:,:,n)] = EBSP_BGCor( ebsp_tile(:,:,n),Settings_Cor );
        end
    else
        ebspBG_tile=ebsp_tile;
    end
    
    clear ebsp_tile

%     %create the BG map array
%     %ebspBG_map=zeros(size(pat_example_bg,1),size(pat_example_bg,2),yn,xn);
%     ebspBG_map=cell(yn,xn);
%     ebsp_map=cell(yn,xn);
%     %ebsp_map=zeros(size(pat_example,1),size(pat_example,2),yn,xn);
% 
%     %reconstruct this into a map
%     %this is not the fastest way to do this, but it is the most readable
%     %it is also only done in RAM

skipped=0;
ebspBG_map=zeros(size(pat_example_bg,1),size(pat_example_bg,2),p_tot);

pmap2=zeros(p_tot,1);

tile.xloaded=zeros(p_tot,1);
tile.yloaded=zeros(p_tot,1);

for xi=1:xn
    for yi=1:yn
        %go from counter to real space
        %go to captured points
        p_i=pmap(yi,xi);

        if isnan(p_i)
            continue
            skipped=skipped+1; %track skipped
        end

        %link to the sorted chunks
        s_n=find(p_sort==p_i);
        
        %put in the map
        
        i=sub2ind([yn,xn],yi,xi)-skipped;
        pmap2(i)=p_i;
        ebspBG_map(:,:,i)=ebspBG_tile(:,:,s_n);
%           ebsp_map{yi,xi}=ebsp_tile(:,:,s_n);

        tile.xloaded(i)=xi;
        tile.yloaded(i)=yi;
    end
end
     clear ebspBG_tile
end

%% Read the EDS data
if PCA_Setup.PCA_EDX==1
    %read the EDX
    % [EDSData1_cor,EDSData1_raw ] = bReadEDX(MapInfo.EBSPData,1,Settings_Cor.channum);

    EDS_cor_block=zeros(Settings_Cor.channum,p_tot);
    EDS_raw_block=zeros(Settings_Cor.channum,p_tot);
    %read the blocks 
    for n=1:numel(block_ends)-1
        [ EDS_cor_block(:,[block_ends(n)+1]:block_ends(n+1)),EDS_raw_block(:,[block_ends(n)+1]:block_ends(n+1))] = bReadEDX(MapInfo.EBSPData,[(p_sort(block_ends(n)+1)) (block_ends(n+1)-block_ends(n))],Settings_Cor.channum);
    end

    %create the BG map array
    EDS_cor_map=zeros(Settings_Cor.channum,p_tot);
    i=1;
    
    %reconstruct this into a map
    %this is not the fastest way to do this, but it is the most readable
    %it is also only done in RAM
    for xi=1:xn
        for yi=1:yn
            %go from counter to real space
            y_r=y_grid(yi,xi);
            x_r=x_grid(yi,xi);
            %go to captured points
            p_i=Data_InputMap.PMap(y_r,x_r);
            
            if isnan(p_i)
                continue
            end

            %link to the sorted chunks
            s_n=find(p_sort==p_i);

            %put in the map
            EDS_cor_map(:,i)=EDS_cor_block(:,s_n);
            
            i=i+1;
        end
    end
    
    clear EDS_cor_block EDS_raw_block
    
end

%% Convolve with the spatial kernel if required

if PCA_Setup.SpatialKernel==1
    
    ebspBG_map=reshape(ebspBG_map,[size(pat_example_bg,1),size(pat_example_bg,2),yn,xn]);
    
    %clear ebspBG_tile %for space
    
    pTime('Applying a spatial kernel',t1);
    
    %Generate the kernela nd get a shorthand for some info
    r=PCA_Setup.KernelRadius;
    extrapts=floor(r/2);
    kernel=GenerateKernel(PCA_Setup.KernelFunction,r);

    %determine whether to load any additional points based on size of AOI
    %kernel has to be odd anyway
    %'X_chunk' corresponds to extra patterns (and their pattern numbers)
    %for padded out AOIs.
    topchunk_cols=[tile.colstart(tiling_no)-extrapts:1:tile.colfin(tiling_no)+extrapts];
    topchunk_rows=[tile.rowstart(tiling_no)-extrapts:1:tile.rowstart(tiling_no)-1];
    topchunk={topchunk_rows,topchunk_cols};
    
    leftchunk_cols=[tile.colstart(tiling_no)-extrapts:1:tile.colstart(tiling_no)-1];
    leftchunk_rows=[tile.rowstart(tiling_no):1:tile.rowfin(tiling_no)];
    leftchunk={leftchunk_rows,leftchunk_cols};
    
    rightchunk_cols=[tile.colfin(tiling_no)+1:1:tile.colfin(tiling_no)+extrapts];
    rightchunk_rows=[tile.rowstart(tiling_no):1:tile.rowfin(tiling_no)];
    rightchunk={rightchunk_rows,rightchunk_cols};
    
    botchunk_cols=[tile.colstart(tiling_no)-extrapts:1:tile.colfin(tiling_no)+extrapts];
    botchunk_rows=[tile.rowfin(tiling_no)+1:1:tile.rowfin(tiling_no)+extrapts];
    botchunk={botchunk_rows,botchunk_cols};
    
    top=Tile_Extend(topchunk,Data_InputMap);
    bot=Tile_Extend(botchunk,Data_InputMap);
    left=Tile_Extend(leftchunk,Data_InputMap);
    right=Tile_Extend(rightchunk,Data_InputMap);
    
    % load and store the extra required patterns
    if PCA_Setup.PCA_EBSD==1
%         top_pats=zeros(Settings_Cor.size,Settings_Cor.size,size(top,1),size(top,2));
%         left_pats=zeros(Settings_Cor.size,Settings_Cor.size,size(left,1),size(left,2));
%         right_pats=zeros(Settings_Cor.size,Settings_Cor.size,size(right,1),size(right,2));
%         bot_pats=zeros(Settings_Cor.size,Settings_Cor.size,size(bot,1),size(bot,2));

        top_pats=zeros(size(pat_example,1),size(pat_example,2),size(top,1),size(top,2));
        left_pats=zeros(size(pat_example,1),size(pat_example,2),size(left,1),size(left,2));
        right_pats=zeros(size(pat_example,1),size(pat_example,2),size(right,1),size(right,2));
        bot_pats=zeros(size(pat_example,1),size(pat_example,2),size(bot,1),size(bot,2));
        
        for i=1:size(top_pats,3)
            for j=1:size(top_pats,4)
                pat=bReadEBSP(MapInfo.EBSPData,top(i,j));
                %refpat_cor=EBSP_BGCor(pat,Settings_Cor);
                top_pats(:,:,i,j)=pat;
            end
        end
        
        for i=1:size(left_pats,3)
            for j=1:size(left_pats,4)
                pat=bReadEBSP(MapInfo.EBSPData,left(i,j));
                %refpat_cor=EBSP_BGCor(pat,Settings_Cor);
                left_pats(:,:,i,j)=pat;
            end
        end
        
        for i=1:size(right_pats,3)
            for j=1:size(right_pats,4)
                pat=bReadEBSP(MapInfo.EBSPData,right(i,j));
                %refpat_cor=EBSP_BGCor(pat,Settings_Cor);
                right_pats(:,:,i,j)=pat;
            end
        end
        
        for i=1:size(bot_pats,3)
            for j=1:size(bot_pats,4)
                pat=bReadEBSP(MapInfo.EBSPData,bot(i,j));
                %refpat_cor=EBSP_BGCor(pat,Settings_Cor);
                bot_pats(:,:,i,j)=pat;
            end
        end
        
        %fill out the AOI with extra patterns
        new_ebsp_array=cat(4,left_pats,ebspBG_map,right_pats);
        new_ebsp_array=cat(3,top_pats,new_ebsp_array,bot_pats);
        
        clear ebspBG_map
        %convolve for each pixel in the EBSP
%         for i=1:size(pat_example,1)
%             parfor j=1:size(pat_example,2)
%                 convolution=conv2(squeeze(new_ebsp_array(i,j,:,:)),kernel,'same');
%                 ebsp_map(i,j,:,:)=convolution(extrapts+1:end-extrapts,extrapts+1:end-extrapts);
%             end
%         end
%         
%         for i=1:size(ebsp_map,3)
%             parfor j=1:size(ebsp_map,4)
%                 ebspBG_map(:,:,i,j)=EBSP_BGCor(ebsp_map(:,:,i,j),Settings_Cor);
%             end
%         end

%       % do spatial averaging with a looped sum
        kernel_x=repmat([1:1:r]-ceil(r/2),[r,1]);
        kernel_y=kernel_x';
        
        kernel_linear=reshape(kernel,[],1);
        kernel_x_lin=reshape(kernel_x,[],1);
        kernel_y_lin=reshape(kernel_y,[],1);
        
        %store the corrected patterns (in slightly the wrong shape +
        %extrapts) in ebspBG_map
        for i=extrapts+1:yn+extrapts
            for j=extrapts+1:xn+extrapts
                centralpat=zeros(size(pat_example,1),size(pat_example,2));
                for k=1:length(kernel_linear)
                    centralpat=centralpat+new_ebsp_array(:,:,i+kernel_x_lin(k),j+kernel_y_lin(k)).*kernel_linear(k);
                end
                ebspBG_map(:,:,i-extrapts,j-extrapts)=EBSP_BGCor(centralpat,Settings_Cor);
            end
        end
        
        
    end
    
    % load the extra required spectra
    if PCA_Setup.PCA_EDX==1
        
        EDS_cor_map=reshape(EDS_cor_map,[Settings_Cor.channum,yn,xn]);
        
        top_spec_cor=zeros(Settings_Cor.channum,size(top,1),size(top,2));
        %top_spec_raw=zeros(Settings_Cor.channum,size(top,1),size(top,2));
        
        left_spec_cor=zeros(Settings_Cor.channum,size(left,1),size(left,2));
        %left_spec_raw=zeros(Settings_Cor.channum,size(left,1),size(left,2));
        
        right_spec_cor=zeros(Settings_Cor.channum,size(right,1),size(right,2));
        %right_spec_raw=zeros(Settings_Cor.channum,size(right,1),size(right,2));
        
        bot_spec_cor=zeros(Settings_Cor.channum,size(bot,1),size(bot,2));
        %bot_spec_raw=zeros(Settings_Cor.channum,size(bot,1),size(bot,2));
        
        for i=1:size(top_spec_cor,2)
            for j=1:size(top_spec_cor,3)
                [spec_cor,spec_raw]=bReadEDX(MapInfo.EBSPData,top(i,j),Settings_Cor.channum);
                top_spec_cor(:,i,j)=spec_cor;
                %top_spec_raw(:,i,j)=spec_raw;
            end
        end
        
        for i=1:size(left_spec_cor,2)
            for j=1:size(left_spec_cor,3)
                [spec_cor,spec_raw]=bReadEDX(MapInfo.EBSPData,left(i,j),Settings_Cor.channum);
                left_spec_cor(:,i,j)=spec_cor;
                %left_spec_raw(:,i,j)=spec_raw;
            end
        end
        
        for i=1:size(right_spec_cor,2)
            for j=1:size(right_spec_cor,3)
                [spec_cor,spec_raw]=bReadEDX(MapInfo.EBSPData,right(i,j),Settings_Cor.channum);
                right_spec_cor(:,i,j)=spec_cor;
                %right_spec_raw(:,i,j)=spec_raw;
            end
        end

        for i=1:size(bot_spec_cor,2)
            for j=1:size(bot_spec_cor,3)
                [spec_cor,spec_raw]=bReadEDX(MapInfo.EBSPData,bot(i,j),Settings_Cor.channum);
                bot_spec_cor(:,i,j)=spec_cor;
                %bot_spec_raw(:,i,j)=spec_raw;
            end
        end
        
    
        % fill out the array with extra required spectra
        new_spec_cor_array=cat(3,left_spec_cor,EDS_cor_map,right_spec_cor);
        %new_spec_raw_array=cat(3,left_spec_raw,EDS_raw_map,right_spec_raw);
        
        new_spec_cor_array=cat(2,top_spec_cor,new_spec_cor_array,bot_spec_cor);
        %new_spec_raw_array=cat(2,top_spec_raw,new_spec_raw_array,bot_spec_raw);
        
        clear EDS_cor_map 
        
        %convolve for each pixel in the spectrum
%         parfor i=1:Settings_Cor.channum
%             convolution_cor=conv2(squeeze(new_spec_cor_array(i,:,:)),kernel,'same');
%             convolution_raw=conv2(squeeze(new_spec_raw_array(i,:,:)),kernel,'same');
%             EDS_cor_map(i,:,:)=convolution_cor(extrapts+1:end-extrapts,extrapts+1:end-extrapts);
%             EDS_raw_map(i,:,:)=convolution_raw(extrapts+1:end-extrapts,extrapts+1:end-extrapts);
%         end

        for i=extrapts+1:yn+extrapts
            for j=extrapts+1:xn+extrapts
                spec_cor=zeros(Settings_Cor.channum,1);
                %spec_raw=zeros(Settings_Cor.channum,1);
                for k=1:length(kernel_linear)
                    spec_cor=spec_cor+new_spec_cor_array(:,i+kernel_x_lin(k),j+kernel_y_lin(k)).*kernel_linear(k);
                    %spec_raw=spec_raw+new_spec_raw_array(:,i+kernel_x_lin(k),j+kernel_y_lin(k)).*kernel_linear(k);
                end
                EDS_cor_map(:,i-extrapts,j-extrapts)=spec_cor;
                %EDS_raw_map(:,i-extrapts,j-extrapts)=spec_raw;
            end
        end
        
    end 
     
end


%% add EBSD data to output
if PCA_Setup.PCA_EBSD==1
    PatternArray=reshape(ebspBG_map,EBSD_Info.PatSizeH*EBSD_Info.PatSizeW,[]);
    Data.Patterns=PatternArray;
    Data.Patterns_reshaped=reshape(Data.Patterns,size(Data.Patterns,1),size(Data.Patterns,2)*size(Data.Patterns,3));
    Data.Patterns_reshaped_norm=Data.Patterns_reshaped./(std(Data.Patterns_reshaped)); %force a normalisation by stdev even if this wasn't in settings_cor
end

%%EDS data normalisation
if PCA_Setup.PCA_EDX==1
EDS_cor_map=reshape(EDS_cor_map(1:Settings_Cor.channum,:),[Settings_Cor.channum,size(EDS_cor_map,2)*size(EDS_cor_map,3)]);
for a=1:size(EDS_cor_map,2)
        integ(a)=sum(EDS_cor_map(:,a)); %normalise by intensity
        EDS_cor_map(:,a)=EDS_cor_map(:,a)./integ(a);

end
EDS_cor_map=EDS_cor_map./std(EDS_cor_map); %normalise by standard deviation.

%add to data output
%Data.EDSData_cor_map=permute(EDS_cor_map,[2,3,1]);
Data.EDSData_cor_normvector=EDS_cor_map;
%Data.EDSData_raw_map=permute(EDS_raw_map,[2,3,1]);
end

end
