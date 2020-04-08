function [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( screensize, LPTsize )
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
%   Detailed explanation goes here

%Fast Fourier transform XCF
%Set up the filters
SettingsXCF.roisize=screensize;
SettingsXCF.roiloc  = [SettingsXCF.roisize/2 SettingsXCF.roisize/2];
SettingsXCF.numroi = 1;
SettingsXCF.filters=[0.5,0.5/2,38,38/2]; %From fiddling these seem to be the best filter settings, but there may be better ones

% SettingsXCF.filters = round([log2(SettingsXCF.roisize)/2,log2(SettingsXCF.roisize)/4,4*log2(SettingsXCF.roisize),2*log2(SettingsXCF.roisize)]);
[SettingsXCF.FFTfilter,SettingsXCF.hfilter] = fFilters(SettingsXCF.roisize,SettingsXCF.filters);

%upsampling mesh - keep small for speed
SettingsXCF.mesh=3;
% SettingsXCF.mesh=1;

%set up the Z rotation LPT settings
[ correction,SettingsXCF2 ] = LPT_Settings( LPTsize,SettingsXCF );

%pack into one output variable
SettingsXCF.roisize2=SettingsXCF2.roisize;
SettingsXCF.roiloc2=SettingsXCF2.roiloc;
SettingsXCF.numroi2=SettingsXCF2.numroi;
SettingsXCF.filters2=SettingsXCF2.filters;
SettingsXCF.FFtfilter2=SettingsXCF2.FFTfilter;
SettingsXCF.hfilter2=SettingsXCF2.hfilter;
SettingsXCF.mesh=SettingsXCF2.mesh;

SettingsXCF.xlim=correction.xlim;
SettingsXCF.xrange=correction.xrange;
SettingsXCF.yroi_num=correction.yroi_num;
SettingsXCF.ycen=correction.ycen;

end

