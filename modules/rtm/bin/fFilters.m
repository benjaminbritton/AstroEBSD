function [FFTfilter,hfilter] = fFilters(roisize,fpassset)
%XCOR_FILTERS Creates filters for EBSD XCorrelation
%   COPYRIGHT TBB 2009
%   INPUTS
%   roisize = roisize, ie 128,256 etc.
%   fpasset = [High Pass Cut Off, High Pass Width, Low Pass Cut Off, Low Pass Width]
%
%   OUTPUTS
%   fftfilter - non idea (gaussian) band pass filter
%   hfilter   - hanning filter
lcutoff=fpassset(3);
lwidth=fpassset(4)/2;
hcutoff=fpassset(1);
hwidth=fpassset(2)/2;

if lcutoff < hcutoff
    error('The low pass filter is smaller than the high pass filter');
end

%generate an x and y grid
u=1:roisize;
[meshv,meshu]=meshgrid(u,u);
meshvf=meshv-roisize/2-0.5;
meshuf=meshu-roisize/2-0.5;

%create the Hfilter
hfilter=cos(pi.*(meshuf)/roisize).*cos(pi.*(meshvf)/roisize);

%create the FFTfilter

distf=sqrt(meshvf.*meshvf+meshuf.*meshuf);

%lowpass
lFFTfilter=exp(-((distf-lcutoff)/(sqrt(2)*lwidth/2)).^2);
lFFTfilter(distf>(lcutoff+2*lwidth))=0;
lFFTfilter(distf<lcutoff)=1;
%highpass
hFFTfilter=exp(-((hcutoff-distf)/(sqrt(2)*hwidth/2)).^2);
hFFTfilter(distf<(hcutoff-2*hwidth))=0;
hFFTfilter(distf>hcutoff)=1;

%combine
FFTfilter=hFFTfilter.*lFFTfilter;
FFTfilter=fftshift(FFTfilter);
end