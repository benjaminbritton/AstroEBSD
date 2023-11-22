%% This is a function that incorporates single hdf file processing for batch processing.

% This script should only be called as a function, and not by itself.

% hdfpatternprocess.m
% This function read and generate hdf5 images, plus primitive image
% processing for DED-captured patterns.
% This is better for instant visualizing or exporting patterns. To blend into
% workflow (indexing), consider using hdfpatternread and then AstroEBSD's 
% EBSP_BGCor function!

% This script uses the most standard Gaussian and gamma filters, which are generally
% suitable for EBSD patterns. Dedicated functions will be available for TKD
% patterns.

% Created by Tianbi Zhang on April 25, 2022
% Last edited by Tianbi on Nov 15, 2023
% Edited by T Ben Britton 
% Requires the MATLAB image processing toolbox!

function [event_data_matrix_cor, bgcor_pattern] = hdfpatternprocess(path, file, gfilts, gamma, outputoption, refoption)
%% This is the main function

%% Toolbox locations for AstroEBSD and MTEX
% Change the following to your AstroEBSD and MTEX directories
% Here we should assume that the user have initialized AstroEBSD in the
% parent script!

%% HDF5 file I/O
h5filename = strcat( file, '.h5');
h5file = fullfile(path, h5filename);

% Read width and height values from the HDF5 file. Both should be 256
% because of the DED geometry
width = h5read(h5file,'/Frame_0/Width');
height = h5read(h5file,'/Frame_0/Height');

% Extract pixel-by-pixel intensity from the "Event" and "iTOT" data.
% Change "Event" to "iToT", "ToT", "ToA" etc. based on data structure.
% Please check the file under HDF5 viewer or use h5info() first.
event_data = double(h5read(h5file,'/Frame_0/SubFrames/Event/Data'));
tot_data= h5read(h5file,'/Frame_0/SubFrames/iToT/Data');

%% Convert data to matrices
% medi_reshape=@(x,w,h) double(reshape(x,w,h));
event_data_matrix=medi_reshape(event_data, width, height);
tot_data_matrix=medi_reshape(tot_data, width, height);

event_data_matrix_cor = imcrop(event_data_matrix,[2 2 251 251]);

% we need to rotate 90 degrees for the EBSD geometry.
if refoption == true
    event_data_matrix_cor = rot90(event_data_matrix_cor,3);
end

%% Image processing
% Define 2x2 Median filter, threshold, dead pixel correction
event_data_matrix_cor = rawhdfcor(event_data_matrix_cor, tot_data_matrix);

%% Background subtraction
% To run this, need to install the statistics and machine learning toolbox!
% From the preamble of EBSP_BGCor
%background correction
Settings_Cor.gfilt=1; %use a low pass filter
Settings_Cor.gfilt_s = gfilts; %low pass filter sigma
%gamma correction
Settings_Cor.gamma = 1;
Settings_Cor.gamma_v = gamma;

% bgcor_pattern = event_data_matrix_cor;
[bgcor_pattern, ~] = EBSP_BGCor(event_data_matrix_cor, Settings_Cor);

bgcor_pattern = normalizeto1(bgcor_pattern);
bgcor_pattern = imadjust(bgcor_pattern,[],[],gamma);

event_data_matrix_cor = normalizeto16bit(event_data_matrix_cor);
bgcor_pattern = normalizeto16bit(bgcor_pattern);

%% Optional - export the patterns

if outputoption == true
bgcor_file_name = strcat(file, '_bgcor.tif');
outputbgcorfile = fullfile(path, bgcor_file_name);
imwrite(uint16(flipud(bgcor_pattern)), outputbgcorfile);

% You can elect to uotput to another image format. Below is an example
% bgcor_file_name_png = strcat(file, '_bgcor.png');
% outputbgcorfile_png = fullfile(path, bgcor_file_name_png);
% imwrite(uint16(flipud(bgcor_pattern)), outputbgcorfile_png);

outputfile = strcat(file, '.tif');
outputfilename = fullfile(path, outputfile);
imwrite(uint16(flipud(event_data_matrix_cor)), outputfilename);
else
    figure;
    subplot(2,1,1);
    imagesc(event_data_matrix_cor); axis xy; axis image; axis off; colormap('gray');
    subplot(2,1,2);
    imagesc(bgcor_pattern); axis xy; axis image; axis off; colormap('gray');
end

end

%% helper functions
function event_matrix = medi_reshape(vector_in, width, height)
event_matrix = double(reshape(vector_in,width,height));
% event_matrix = event_matrix';
end

function corrected_event_matrix = rawhdfcor(event_matrix_in, tot_matrix_in)
event_data_matrix_filter = medfilt2(event_matrix_in, [3 3]); % Define 2x2 Median filter
med_dif = event_matrix_in - event_data_matrix_filter; % compare the original and filtered data sets
med_thresh_pos =  6*std(med_dif(:)); % define a threshold difference
med_thresh_neg = -med_thresh_pos;
% 
% 
% % apply the median filter
corrected_event_matrix = event_matrix_in;
corrected_event_matrix(event_matrix_in == 0) = event_data_matrix_filter(event_matrix_in == 0);
% corrected_event_matrix(tot_matrix_in == 0) = event_data_matrix_filter(tot_matrix_in == 0);

% For those higher than the upper threshold, or lower than the lower
% threshold, replace the pixel data by the filtered pixel value
corrected_event_matrix(med_dif > med_thresh_pos) = event_data_matrix_filter(med_dif > med_thresh_pos);
corrected_event_matrix(med_dif < med_thresh_neg) = event_data_matrix_filter(med_dif < med_thresh_neg);
% Check the dead pixels based on iTOT data


end

function normalized_matrix = normalizeto1(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
end

function normalized_matrix = normalizeto16bit(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
normalized_matrix = normalized_matrix * (2^16 - 1); % Convert to 16 bit
end
