function [FFTData,data_fill]=fROIEx2(Image1,SettingsXCF)
%extract FFT from ROI
% This code is copyright TBB 2010
% b.britton@imperial.ac.uk

%put together the bounding window sizes
roi_imgx=zeros(SettingsXCF.roisize,SettingsXCF.numroi);
roi_imgy=zeros(SettingsXCF.roisize,SettingsXCF.numroi);


for n=1:SettingsXCF.numroi
    roi_imgy(:,n)=(1:SettingsXCF.roisize)+SettingsXCF.roiloc(n,2)-SettingsXCF.roisize/2;
    roi_imgx(:,n)=(1:SettingsXCF.roisize)+SettingsXCF.roiloc(n,1)-SettingsXCF.roisize/2;
end

data_fill=[1:(SettingsXCF.filters(3)+SettingsXCF.filters(4)),SettingsXCF.roisize-(SettingsXCF.filters(3)+SettingsXCF.filters(4)-1):SettingsXCF.roisize];
data_fillsize=length(data_fill);
SettingsXCF.FFTfilter=SettingsXCF.FFTfilter(data_fill,data_fill);

FFTData=complex(zeros(data_fillsize,data_fillsize,SettingsXCF.numroi),zeros(data_fillsize,data_fillsize,SettingsXCF.numroi));

for n=1:SettingsXCF.numroi
    roi_ref_temp=Image1(roi_imgy(:,n),roi_imgx(:,n));
    roi_ref_temp=SettingsXCF.hfilter.*(roi_ref_temp-nanmean(roi_ref_temp(:)))/nanstd(roi_ref_temp(:));
    
    roi_ref_tempf=fft2(roi_ref_temp);
    %reduce the data range
    FFTData(:,:,n)=SettingsXCF.FFTfilter.*roi_ref_tempf(data_fill,data_fill);
end


end