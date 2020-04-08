% Tile dimensions and locations - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [MapInfo]=TilingSetup(MapInfo,crop_factor)

    %Define cropped dimensions
    MapInfo.cropped_width=floor(MapInfo.MapSize2./crop_factor);
    MapInfo.cropped_height=floor(MapInfo.MapSize1./crop_factor);
    
    %Adjust cropping window so they are all the same shape
    a=1;
    while a==1
        if MapInfo.cropped_width.*crop_factor<=MapInfo.MapSize2 && MapInfo.cropped_height.*crop_factor<=MapInfo.MapSize1
            break
        end
        
        if MapInfo.cropped_width.*crop_factor>MapInfo.MapSize2
            MapInfo.cropped_width=MapInfo.cropped_width-1;
        end
        
        if MapInfo.cropped_height.*crop_factor>MapInfo.MapSize1
            MapInfo.cropped_height=MapInfo.cropped_height-1;
        end
    end
end