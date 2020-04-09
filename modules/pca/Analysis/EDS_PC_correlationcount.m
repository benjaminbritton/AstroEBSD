%PC EDS spectra correlation counts - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Quality matching map of EDS spectra
%Produces correlation quality with all other PCs (in quadrature), count
%skipped if too high a match, and count non-skipped.

total_ret_comps=sum(tile.ret_comps(1:end));
qual=zeros(total_ret_comps,total_ret_comps);

for a=1:total_ret_comps
spec1=tile.spec_reshaped(:,a);
    
quality(a)=0;
skipped(a)=0;
count(a)=0;
for b=1:total_ret_comps   
    spec2=tile.spec_reshaped(:,b);
    qual(a,b)=max(xcorr(spec1,spec2,'coeff'));%./test;
    
    if qual(a,b)>0.98
        skipped(a)=skipped(a)+1;
        count(a)=count(a)+0;
        continue
    end
    
    quality(a)=(quality(a)+qual(a,b).^2);
    count(a)=count(a)+1;
    %EDScorrqual(a,b)=quality;
end
disp(a)
end
EDSWindowing.qual=qual;
EDSWindowing.count=count;
EDSWindowing.skipped=skipped;
EDSWindowing.EDScorrqual=quality;
clear quality count qual skipped spec1 spec2
%EDSWindowing.EDScorrqual2=EDSWindowing.EDScorrqual.*abs(1-eye(size(EDSWindowing.EDScorrqual))); %set all diagonal elements to zero
%EDSWindowing.corrqualsum=sqrt(sum(EDSWindowing.EDScorrqual.^2));