%  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function []=ElementMaps_combined(tile,quant,elements,comb,printing)

cmap=cbrewer('qual','Set1',length(elements));


for el=1:size(comb,1)
    
    r=comb(el,1);
    g=comb(el,2);
    b=comb(el,3);
    
mapr=zeros(size(tile.map_reshaped)); 
mapg=zeros(size(tile.map_reshaped)); 
mapb=zeros(size(tile.map_reshaped)); 
for pc=1:size(quant,2) 
mapr((tile.map_reshaped==pc))=quant(r,pc);
mapg((tile.map_reshaped==pc))=quant(g,pc);
mapb((tile.map_reshaped==pc))=quant(b,pc);
end

mapr=(uint8(255.*(mapr./max(max(mapr)))));
mapg=(uint8(255.*mapg./max(max(mapg))));
mapb=(uint8(255.*mapb./max(max(mapb))));

maprgb(:,:,1)=mapr;
maprgb(:,:,2)=mapg;
maprgb(:,:,3)=mapb;

imshow(maprgb)

if printing==1
savefig([elements{r},elements{g},elements{b}])
print(gcf,[elements{r},elements{g},elements{b}],'-dpng','-r500');
end

end
    
end