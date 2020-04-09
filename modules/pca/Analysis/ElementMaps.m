% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function []=ElementMaps(tile,quant,elements,printing)

cmap=cbrewer('qual','Set1',length(elements));


for el=1:size(quant,1)
    
map=zeros(size(tile.map_reshaped)); 
for pc=1:size(quant,2)
map((tile.map_reshaped==pc))=quant(el,pc);
end

map=map./max(max(map));
imagesc(map)
colormap( cmap(el,:).* [0:0.01:1]' );

if printing==1
savefig([elements{el}])
print(gcf,[elements{el}],'-dpng','-r500');
end

end
    
end