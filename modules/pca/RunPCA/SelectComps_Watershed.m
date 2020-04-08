% Watershed component selection - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%This script uses a series of filters to label a radon quality map using
%the watershed transform.

%Grab the radon map
pq=Plots.RadonCropped(:,:,tiling_no);
%Guided filter
pqe=imguidedfilter(pq,'NeighborhoodSize',[10,10],'DegreeOfSmoothing',0.001*diff(getrangefromclass(pq)).^2);
imagesc(pqe)
%%
%STD filter then convolve with a blurring kernel
pq2=conv2(stdfilt(pqe),[1,1;1,1]);
%Take the complement
pqfilt=imcomplement(255.*pq2);
%Gaussian filter
pqfilt2=imgaussfilt(pqfilt);
I2 = imcomplement(pqfilt2);
%Do watershed transform
I3 = imhmin(I2,0.33333*mean(mean(std(I2)))); 
L = watershed(I3);
imagesc(L)
Watershed_selcomps=(max(max(sort(L))));

% extra_comps=floor(0.2.*ret_comps);
% if extra_comps > 5
%     extra_comps=5;
% end
% ret_comps=ret_comps+extra_comps;

tile.Watershed_selcomps(tiling_no)=Watershed_selcomps;
tile.RadonWatershed(:,:,tiling_no)=L;
clear pq pq2 pqfilt pqfilt2 I2 I3 L pqe extra_comps