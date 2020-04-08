function [FFTData,data_fill]=fROIEx(Image1,XCF_hfilter,XCF_FFTfilter,XCF_roiloc,XCF_filters,XCF_numroi,XCF_roisize)

%extract FFT from ROI

% This code is copyright TBB 2010
% do not share without express permission.
% b.britton@imperial.ac.uk

%put together the bounding window sizes
roi_imgx=zeros(XCF_roisize,XCF_numroi);
roi_imgy=zeros(XCF_roisize,XCF_numroi);


for n=1:XCF_numroi
    roi_imgy(:,n)=(1:XCF_roisize)+XCF_roiloc(n,2)-XCF_roisize/2;
    roi_imgx(:,n)=(1:XCF_roisize)+XCF_roiloc(n,1)-XCF_roisize/2;
end

data_fill=[1:(XCF_filters(3)+XCF_filters(4)),XCF_roisize-(XCF_filters(3)+XCF_filters(4)-1):XCF_roisize];
data_fillsize=length(data_fill);
XCF_FFTfilter=XCF_FFTfilter(data_fill,data_fill);

FFTData=complex(zeros(data_fillsize,data_fillsize,XCF_numroi),zeros(data_fillsize,data_fillsize,XCF_numroi));

for n=1:XCF_numroi
    roi_ref_temp=Image1(roi_imgy(:,n),roi_imgx(:,n));
    roi_ref_temp=XCF_hfilter.*(roi_ref_temp-nanmean(roi_ref_temp(:)))/nanstd(roi_ref_temp(:));
    
    roi_ref_tempf=fft2(roi_ref_temp);
    %reduce the data range
    FFTData(:,:,n)=XCF_FFTfilter.*roi_ref_tempf(data_fill,data_fill);
end


end