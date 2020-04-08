function [ correction,SettingsXCF2 ] = LPT_Settings( LPTsize,SettingsXCF )

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
correction.xlim=150:406;
correction.xlim=round(LPTsize*correction.xlim/500);
correction.xrange=range(correction.xlim)+1; %Bits for the smaller ROIs
% FFT reference

SettingsXCF2=SettingsXCF;
SettingsXCF2.mesh = 3;
SettingsXCF2.roisize=correction.xrange; %Settings for the FFT in the z correction, these need to be different to the ones for the FFT indexing
% 
%set up the filters
SettingsXCF.filters=[0.5,0.5/2,38,38/2]; %From fiddling these seem to be the best filter settings, but there may be better ones
% SettingsXCF2.filters = round([log2(SettingsXCF2.roisize)/2,log2(SettingsXCF2.roisize)/4,4*log2(SettingsXCF2.roisize),2*log2(SettingsXCF2.roisize)]); %playing with these filters might help thinking of it
[SettingsXCF2.FFTfilter,SettingsXCF2.hfilter] = fFilters(SettingsXCF2.roisize,SettingsXCF2.filters);

correction.yroi_num=3;
correction.ycen=round(linspace(1,LPTsize-correction.xrange,correction.yroi_num));

end

