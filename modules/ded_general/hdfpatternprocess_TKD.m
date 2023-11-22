%% This is a function that incorporates single hdf file processing for batch processing.

% This script should only be called as a function, and not by itself.

% hdfpatternprocess_satcor.m
% This function read and display hdf5 images, plus primitive image
% processing for DED-captured patterns. Workflow is as the following:
% (1) The scriptreads in data in the hdf file
% (2) The script identifies the saturation region based on pixel intensity
% (user input) 
% (3) A dead pixel correction step using iToT data and a 3x3 median filter
% (4) A 25x25 median filter of the raw data is used to flatten the raw
% pattern
% (5) Saturated region on the flattened pattern is replaced by random noise
% with mean and stdev the same as the rest of the flattened pattern.
% (6) An optional gamma filter is available. If you are not sure, leave as
% gamma = 1.
% (7) You can opt to have images written to files, or just display the raw and
% background corrected pattern on a MATLAB window through 'outputoption'.

% The saturation index must be carefully chosen so useful features like
% diffraction spots are not unintentionally wiped out!!

% For off-axis TKD or EBSD patterns, please consider using the primitive version
% "hdfpatternprocess" which uses a simpler Gaussian filter without saturation
% correction.
% Created by Tianbi Zhang on Aug 22, 2023
% Edited by T Ben Britton 
% Requires the MATLAB image processing toolbox!

function [event_data_matrix_cor, bgcor_pattern] = hdfpatternprocess_satcor(path, file, gamma, satcor_satindex, outputoption)
%% HDF5 file I/O
rawfile = strcat( path, file);
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
medi_reshape=@(x,w,h) double(reshape(x,w,h));
event_data_matrix=medi_reshape(event_data, width, height)';
tot_data_matrix=medi_reshape(tot_data, width, height)';

event_data_matrix = imcrop(event_data_matrix,[3 3 251 251]);
tot_data_matrix = imcrop(tot_data_matrix,[3 3 251 251]);

event_data_matrix_cor = rawhdfcor(event_data_matrix, tot_data_matrix, satcor_satindex);

new_event_data = reshape(event_data_matrix,[],1);
satregion_ind = find(new_event_data >= satcor_satindex);
unsatregion_ind = find(new_event_data < satcor_satindex);

%% Flatten the pattern
medf = medfilt2(event_data_matrix_cor,[25 25],'symmetric');
bgcor_pattern = event_data_matrix_cor ./ medf;

cor_vec = reshape(bgcor_pattern, [],1);

%% Create and fill in the noise
pat_mean = mean(cor_vec(unsatregion_ind));
pat_stdev = std(cor_vec(unsatregion_ind));
y = pat_stdev.*randn(size(satregion_ind)) + pat_mean;

bgcor_pattern(satregion_ind) = y;

bgcor_pattern = normalizeto1(bgcor_pattern);
bgcor_pattern = imadjust(bgcor_pattern,[],[],gamma);
bgcor_pattern = normalizeto16bit(bgcor_pattern);
event_data_matrix_cor = normalizeto16bit(event_data_matrix_cor);

%% Plot background corrected pattern side by side with raw pattern
if outputoption == false
figure
title(rawfile);
hold on;
subplot(1,2,1);
imagesc(event_data_matrix_cor);
colormap('gray');
axis image;
axis xy;
title("Raw Event");
subplot(1,2,2);
imagesc(bgcor_pattern);
hold on;
colormap('gray');
axis image;
axis xy;
title("Raw + bgcor");
pause(1);
else
bgcor_file_name = strcat(file, '_bgcor.tif');
outputbgcorfile = fullfile(path, bgcor_file_name);
imwrite(uint16(flipud(bgcor_pattern)), outputbgcorfile);

% bgcor_file_name_png = strcat(file, '_bgcor.png');
% outputbgcorfile_png = fullfile(path, bgcor_file_name_png);
% imwrite(uint16(flipud(bgcor_pattern)), outputbgcorfile_png);

outputfile = strcat(file, '.tif');
outputfilename = fullfile(path, outputfile);
imwrite(uint16(flipud(event_data_matrix_cor)), outputfilename);
end

end

%% helper functions

function corrected_event_matrix = rawhdfcor(event_matrix_in, tot_matrix_in, satcor_satindex)
corrected_event_matrix = event_matrix_in;
event_data_matrix_filter = medfilt2(event_matrix_in, [2 2]); % Define 2x2 Median filter
corrected_event_matrix(event_matrix_in == 0) = event_data_matrix_filter(event_matrix_in == 0);

med_dif = event_matrix_in - event_data_matrix_filter; % compare the original and filtered data sets
med_thresh_pos = satcor_satindex/10; % define a threshold difference
med_thresh_neg = -med_thresh_pos;

corrected_event_matrix(med_dif > med_thresh_pos) = event_data_matrix_filter(med_dif > med_thresh_pos);
corrected_event_matrix(med_dif < med_thresh_neg) = event_data_matrix_filter(med_dif < med_thresh_neg);
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
