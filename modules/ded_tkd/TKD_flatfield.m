function [fused_flat,fusedpattern,normfactor, ref_curve] = TKD_flatfield(rawdata,normalizeddata,angle_matrix,fillingdomains,numfiles)
%TKD_FLATFIELD fit to the radial function, and then flat field with it


[normfactor, ref_curve] = TKD_referencecurvefit(rawdata,normalizeddata,angle_matrix,numfiles);

% Stitch the pattern
[fusedpattern] = TKD_patternfuse(rawdata,normalizeddata,fillingdomains);

%% Flat field
npix=sqrt(size(rawdata,1));

normfactor = reshape(normfactor, [npix npix]);
fused_flat = fusedpattern ./ normfactor;

end

