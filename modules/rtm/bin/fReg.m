 function [output] = fReg( buf1ft,buf2ft,XCF_roisize,XCF_mesh,data_fill)
%FREG Register two FFTs to subpixel accuracy
%a reduced form of the code submitted to the matlab file exchange
%http://www.mathworks.co.uk/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation
%
%reported in the literature:
%Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup
%"Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
%
%also modded to handle the filtered FFT sizes
%TBB 2012

CC=zeros(XCF_roisize*2);

red_roisize=length(data_fill);

CC(XCF_roisize+1-fix(red_roisize/2):XCF_roisize+1+fix((red_roisize-1)/2),XCF_roisize+1-fix(red_roisize/2):XCF_roisize+1+fix((red_roisize-1)/2)) = ...
    fftshift((buf1ft).*conj(buf2ft));

% Compute crosscorrelation and locate the peak
CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
[max1,loc1] = max(CC);
[~,loc2] = max(max1);
rloc=loc1(loc2);cloc=loc2;

% Obtain shift in original pixel grid from the position of the crosscorrelation peak
XCF_roisize2=2*XCF_roisize;

if rloc > XCF_roisize
    row_shift = rloc - XCF_roisize2 - 1;
else
    row_shift = rloc - 1;
end
if cloc > XCF_roisize
    col_shift = cloc - XCF_roisize2 - 1;
else
    col_shift = cloc - 1;
end
row_shift=row_shift/2;
col_shift=col_shift/2;

%%% DFT computation %%%
% Initial shift estimate in upsampled grid
row_shift = round(row_shift*XCF_mesh)/XCF_mesh;
col_shift = round(col_shift*XCF_mesh)/XCF_mesh;
dftshift = fix(ceil(XCF_mesh*1.5)/2); %% Center of output array at dftshift+1

% Matrix multiply DFT around the current shift estimate
roff=dftshift-row_shift*XCF_mesh;
coff=dftshift-col_shift*XCF_mesh;

% Compute kernels and obtain DFT by matrix products
prefac=-1i*2*pi/(XCF_roisize*XCF_mesh);

%speed up kernel generation for reduced filtered FFT
c_i=ifftshift([0:XCF_roisize-1]).' - floor(XCF_roisize/2) ;
c_i=c_i(data_fill);

r_i=ifftshift([0:XCF_roisize-1]) - floor(XCF_roisize/2);
r_i=r_i(data_fill);

kernc=exp(prefac*( c_i )*( [0:ceil(XCF_mesh*1.5)-1] - coff ));
kernr=exp(prefac*( [0:ceil(XCF_mesh*1.5)-1].' - roff )*( r_i ));
kern=buf2ft.*conj(buf1ft);
CC2=conj(kernr*kern*kernc);

% Locate maximum and map back to original pixel grid
[max1,loc1] = max(CC2);
[CCmax,loc2] = max(max1);
CCmax=abs(CCmax);

rloc = loc1(loc2) - dftshift - 1;
cloc = loc2 - dftshift - 1;
%     CCmax = CC(rloc,cloc);

row_shift = row_shift + rloc/XCF_mesh;
col_shift = col_shift + cloc/XCF_mesh;

BF1=sum(buf1ft(:).*conj(buf1ft(:)));
BF2=sum(buf2ft(:).*conj(buf2ft(:)));
XCF_ph=CCmax./sqrt(BF1*BF2);


output=[col_shift,row_shift,CCmax/(XCF_roisize^2),XCF_ph];

end