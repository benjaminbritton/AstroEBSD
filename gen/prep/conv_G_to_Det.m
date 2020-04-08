function [G_output] = conv_G_to_Det(G_input,Detector_tilt)
%CONV_G_TO_DET Converts a G matrix from the sample coordinate system into the detector using Detector tilt

G_output=zeros(3,3,size(G_input,3));
for n=1:size(G_input,3)
    G_output(:,:,n)=G_input(:,:,n)'*Detector_tilt;
end

