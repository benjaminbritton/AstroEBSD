function [ebsp_tile,ebspBG_tile] = bBlockReadEBSP(xpts,ypts,MapInfo,Settings_Cor)
%BBLOCKREADEBSP Read a block of EBSPs as quick as we can
%background correct if needed

%convert to p_numbers
p_tot=numel(xpts);

%read the 1st pattern to get dimensions
[ pat_example ] = bReadEBSP(MapInfo.EBSPData,MapInfo.Data_InputMap.PMap(ypts(1),xpts(1)));
pwid=size(pat_example,2);
phigh=size(pat_example,1);

%create a list of pattern numbers that correspond to
%the grid positions you want
p_nums=zeros(p_tot,1);
for p=1:p_tot
    p_nums(p)=MapInfo.Data_InputMap.PMap(ypts(p),xpts(p));
end

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

ebsp_tile=zeros(phigh,pwid,p_tot);


% %failed attempt at making things parallel - turns out to be slower,
% likley due to the cell array building
% num_chunk=numel(block_ends)-1;

% ebsp_tile_c=cell(num_chunk,1);
% %read the blocks
% EBSP_data=MapInfo.EBSPData;

% parfor n=1:num_chunk
%     ebsp_tile_c{n}=bReadEBSP(EBSP_data,[(p_sort(block_ends(n)+1)) (block_ends(n+1)-block_ends(n))]);
% end

% for n=1:num_chunk
%     ebsp_tile(:,:,[block_ends(n)+1]:block_ends(n+1))=ebsp_tile_c{n};
% end


% %read the blocks
for n=1:numel(block_ends)-1
    [ ebsp_tile(:,:,[block_ends(n)+1]:block_ends(n+1))] = bReadEBSP(MapInfo.EBSPData,[(p_sort(block_ends(n)+1)) (block_ends(n+1)-block_ends(n))]);
end

%reorder
% ebsp_tile=zeros(size(ebsp_tilet));
iv=zeros(p_tot,1);

for p=1:p_tot
    iv(p)=find(p_sort==p_nums(p));
end

%swap the order of this data
ebsp_tile=ebsp_tile(:,:,iv);

if nargout > 1 %BG correct if requested
    [ pat_example_bg] = EBSP_BGCor( pat_example,Settings_Cor );
    ebspBG_tile=zeros(size(pat_example_bg,1),size(pat_example_bg,2));
  
    parfor p=1:p_tot
        [ ebspBG_tile(:,:,p)] = EBSP_BGCor( ebsp_tile(:,:,p),Settings_Cor );
    end
end

end

