function [ Best_G, x_shift, y_shift, peak_height ] = index_FFT2( SettingsXCF, P_ref, G_test,FFTData_test_all )
%INDEX_FFT Summary of this function goes here

% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

library_size=size(G_test,3);

%FFT the refernece image
[FFTData_Ref,XCF_data_fill]  =fROIEx(P_ref,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize);

%Fancy Parallel computig shiz

%do the correaltion
RegOut1=zeros(library_size,4);

for n=1:library_size
    %Cross correlate and find the peak
    [RegOut1(n,:)] = fReg( FFTData_Ref,FFTData_test_all(:,:,n),SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
end

 [peak_height,max_index]=max(RegOut1(:,4));

Best_G = G_test(:,:, max_index);
x_shift = RegOut1(max_index, 1);
y_shift = RegOut1(max_index, 2);
end

