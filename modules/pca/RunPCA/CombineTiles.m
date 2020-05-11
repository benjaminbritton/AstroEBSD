% CombineTiles  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Put tiles together
% Generates maps for full AOI based on PCA/VM on each tile
% Stored in tile structure

n=PCA_Setup.crop_factor;

for a=1:n.^2
    
if PCA_Setup.MedFilt==1
tile.map(:,:,a)=medfilt2(tile.map(:,:,a),'symmetric');
end

end

tile.map_reshaped=[];
tile.map_reshaped=tile.map(:,:,1);
tile.map_reshaped2=tile.map(:,:,1);
%tile.watershed_reshaped=tile.RadonWatershed(:,:,1);
tile.loadingmap_reshaped=tile.loadingmap(:,:,1);

if n>1
for iii=2:n.^2
    tile.map_reshaped=[tile.map_reshaped;tile.map(:,:,iii)+sum(tile.ret_comps(1:iii-1))];
    tile.map_reshaped2=[tile.map_reshaped2;tile.map(:,:,iii)];
    %tile.watershed_reshaped=[tile.watershed_reshaped;tile.RadonWatershed(:,:,iii)];
    tile.loadingmap_reshaped=[tile.loadingmap_reshaped;tile.loadingmap(:,:,iii)];
end

cutpoint=MapInfo.cropped_height.*n;

map_reshaped=tile.map_reshaped(1:cutpoint,:);
map_reshaped2=tile.map_reshaped2(1:cutpoint,:);
%watershed_reshaped=tile.watershed_reshaped(1:cutpoint,:);
loadingmap_reshaped=tile.map_reshaped(1:cutpoint,:);

for iii=1:n-1
    map_reshaped=[map_reshaped,tile.map_reshaped(iii.*cutpoint+1:cutpoint.*(iii+1),:)];
    map_reshaped2=[map_reshaped2,tile.map_reshaped2(iii.*cutpoint+1:cutpoint.*(iii+1),:)];
    %watershed_reshaped=[watershed_reshaped,tile.watershed_reshaped(iii.*cutpoint+1:cutpoint.*(iii+1),:)];
    loadingmap_reshaped=[loadingmap_reshaped,tile.loadingmap_reshaped(iii.*cutpoint+1:cutpoint.*(iii+1),:)];
end


tile.map_reshaped=map_reshaped;
tile.map_reshaped2=map_reshaped2;

%ensure anything that was zero before is zero after
tile.map_reshaped(tile.map_reshaped2==0)=0;

%tile.watershed_reshaped=watershed_reshaped;
tile.loadingmap_reshaped=loadingmap_reshaped;
clear map_reshaped cutpoint map_reshaped2 watershed_reshaped loadingmap_reshaped
end

% tile.map_reshaped=[tile.map(:,:,1),sum(tile.ret_comps(1:3))+tile.map(:,:,4),sum(tile.ret_comps(1:6))+tile.map(:,:,7);...
%     sum(tile.ret_comps(1:1))+tile.map(:,:,2),sum(tile.ret_comps(1:4))+tile.map(:,:,5),sum(tile.ret_comps(1:7))+tile.map(:,:,8);...
%     sum(tile.ret_comps(1:2))+tile.map(:,:,3),sum(tile.ret_comps(1:5))+tile.map(:,:,6),sum(tile.ret_comps(1:8))+tile.map(:,:,9)];
% 
% tile.map_reshaped2=[tile.map(:,:,1),tile.map(:,:,4),tile.map(:,:,7);...
%     tile.map(:,:,2),tile.map(:,:,5),tile.map(:,:,8);...
%     tile.map(:,:,3),tile.map(:,:,6),tile.map(:,:,9)];
% 
% tile.watershed_reshaped=[tile.RadonWatershed(:,:,1),tile.RadonWatershed(:,:,4),tile.RadonWatershed(:,:,7);...
%     tile.RadonWatershed(:,:,2),tile.RadonWatershed(:,:,5),tile.RadonWatershed(:,:,8);...
%     tile.RadonWatershed(:,:,3),tile.RadonWatershed(:,:,6),tile.RadonWatershed(:,:,9)];

if PCA_Setup.PCA_EDX==1
    tile.spec_reshaped=reshape(tile.spec,Settings_Cor.channum,sum(tile.ret_comps(1:end)));
end

if PCA_Setup.PCA_EBSD==1
    tile.pat_reshaped=reshape(tile.pat,EBSD_Info.PatSizeH*EBSD_Info.PatSizeW,sum(tile.ret_comps(1:end)));
end

clear iii n