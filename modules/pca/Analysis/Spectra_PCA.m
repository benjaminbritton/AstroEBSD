% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [data,elements]=Spectra_PCA(loc)
cd(loc)

[data,headers]=xlsread('Results');
data=data(1:end,:);
elements=headers(1,2:end);
%spec_list=headers(2:end-3,1);
spec_list=headers(1:end,1);

%search for and delete the 'sum' column
sumind=find(contains(elements,'Sum'));
elements(sumind)=[];
data(:,sumind)=[];

%%
for b=1:length(spec_list)
    speccheck=spec_list{b};
    if speccheck(end-3:end)=='.spx'
        
        spec_order=spec_list{b};
        spec_order=spec_order(10:end-4);
        order(b)=str2num(spec_order);
    end
end

%
[~,idx]=sort(order);
spec_list = spec_list(idx,:);
data=data(idx,:);

end

