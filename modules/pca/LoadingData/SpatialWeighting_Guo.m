% Spatial kernel for PCA  - T P McAuliffe 20/02/2020 for AstroEBSD v2

% This function is unused in the main PCA script - it generates a new testarray using a previous one and a given kernel.
% For our purposes we now apply the kernel upon data loading to the raw EBSPs

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% This is implemented from Spatially Weighted Principal Component Analysis
% for Imaging Classification - Ruixin Guo, Mihye Ahn, Hongtu Zhu - doi:10.1080/10618600.2014.912135

function [tA_h]=SpatialWeighting_Guo(testArray,MapInfo,PCA_Setup)

%Conventionalise testarray
testArray=testArray';

r = PCA_Setup.KernelRadius;

tA_h=zeros(size(testArray));
tA_size=size(testArray);

d=zeros(1,MapInfo.cropped_width.*MapInfo.cropped_height);
%%
parfor i=1:MapInfo.cropped_width.*MapInfo.cropped_height
        w=[];
        d=[];
        %euclids=[];
        
        for j=1:MapInfo.cropped_width.*MapInfo.cropped_height;
            [y1,x1]=ind2sub([MapInfo.cropped_height,MapInfo.cropped_width],i);
            [y2,x2]=ind2sub([MapInfo.cropped_height,MapInfo.cropped_width],j);
            d(j)=sqrt((x2-x1).^2+(y2-y1).^2);
        end
        
        kernel=find(d<=r);
        
        for h=1:length(kernel)
            %euclids(h)=norm(testArray(i,:)-testArray(kernel(h),:));
            w(h) = (1 - (d(kernel(h))./r).^2).^2;%.* exp(-euclids(h)./e_bw);
        end
        
        %normalise w and euclids - now each row of w sums to 1
        w=w./sum(w);
        w_kernel=transpose(repmat(w,tA_size(2),1));
        
        %get local average pattern vectors and insert into data matrix
        x_h=testArray(kernel,:);
        x_ih=sum(x_h.*w_kernel);
        %x_ih=x_ih./length(kernel);
        
        tA_h(i,:)=x_ih;
end

%bring tA_h back to our convention
tA_h=tA_h';

end
