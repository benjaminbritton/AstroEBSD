function [ rotmat_z_correctionFFT ] = Angle_ShiftZ( EBSP_logref,EBSP_loghalfcorrFFT,correction,LPTsize,SettingsXCF2,Rz )
%ANGLE_SHIFTZ Use the LPT to measure the Z rotation for two EBSPs

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

RegOut3=zeros(correction.yroi_num,4);

for n=1:correction.yroi_num
    
    %set up the smaller ROI from the LPT
    correction.yrange=correction.ycen(n)+(0:correction.xrange-1);
    
    %extract the ROI from both images
    correction.ROI_ref=EBSP_logref(correction.yrange,correction.xlim);
    correction.ROI_test=EBSP_loghalfcorrFFT(correction.yrange,correction.xlim);
    
    %FFT these ROIs
    [FFTData_LPT_ref,XCF_LPT_data_fill]  =fROIEx(correction.ROI_ref,SettingsXCF2.hfilter,SettingsXCF2.FFTfilter,[correction.xrange/2 correction.xrange/2],SettingsXCF2.filters,SettingsXCF2.numroi,SettingsXCF2.roisize);
    [FFTData_LPT_test,~]  =fROIEx(correction.ROI_test,SettingsXCF2.hfilter,SettingsXCF2.FFTfilter,[correction.xrange/2 correction.xrange/2],SettingsXCF2.filters,SettingsXCF2.numroi,SettingsXCF2.roisize);
    
    %XCF this ROI
    [RegOut3(n,:)] = fReg( FFTData_LPT_ref,FFTData_LPT_test,SettingsXCF2.roisize,SettingsXCF2.mesh,XCF_LPT_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    
end

%calculate the z rotation
shift.zFFT = mean(RegOut3(:,2));
sf = LPTsize/360;
angle.zcorrFFT = (shift.zFFT/sf)*pi/180;

%create a matrix
rotmat_z_correctionFFT = Rz(angle.zcorrFFT);

end

