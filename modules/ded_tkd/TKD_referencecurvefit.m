function [normfactor, ref_curve] = TKD_referencecurvefit(rawdata,normalizeddata,angle_matrix,numfiles)
%fit to the angular data to get the scattering profile

%% Fit reference curve to a smoothing function
notdeadpixel = find(rawdata(:,1) ~=0);

ref_curve = normalizeddata(:,numfiles);

normfactor = smooth(angle_matrix(:), ref_curve,200,'sgolay',2);

end