%% This is a function that reads a h5 pattern into a matrix variable

% This script should only be called as a function, and not by itself.

% No normalization is done here!!!

% hdfpatternread
% This function reads the event data in a h5 file, reshape it to a matrix 
% and crop it (the output).
% No background correction or image export is done here.
% 
% Created by Tianbi Zhang on April 25, 2022
% Last edited on November 15, 2023
% Edited by T Ben Britton 
% Requires the MATLAB image processing toolbox!

function corrected_rawpattern = hdfpatternread(path, file, refoption)
%% This is the main function

%% Toolbox locations for AstroEBSD and MTEX
% Change the following to your AstroEBSD and MTEX directories
% Here we should assume that the user have initialized AstroEBSD in the
% parent script!

%% HDF5 file I/O
% rawfile = strcat( path, file);
h5filename = strcat( file, '.h5');
h5file = fullfile(path, h5filename);

% Read width and height values from the HDF5 file. Both should be 256
% because of the DED geometry
width = h5read(h5file,'/Frame_0/Width');
height = h5read(h5file,'/Frame_0/Height');

% Extract pixel-by-pixel intensity from the "Event" and "iTOT" data.
% Change "Event" to "iToT", "ToT", "ToA" etc. based on data structure.
% Please check the file under HDF5 viewer or use h5info() first.
event_data = h5read(h5file,'/Frame_0/SubFrames/Event/Data');
tot_data= h5read(h5file,'/Frame_0/SubFrames/iToT/Data');

%% Convert data to matrices
% medi_reshape=@(x,w,h) double(reshape(x,w,h));
rawpattern=medi_reshape(event_data, width, height);
tot_data_matrix=medi_reshape(tot_data, width, height);

corrected_rawpattern = rawhdfcor(rawpattern, tot_data_matrix);

% we need to rotate 90 degrees for the EBSD geometry.
if refoption == true
    corrected_rawpattern = rot90(corrected_rawpattern,3);
end

corrected_rawpattern = imcrop(corrected_rawpattern,[2 2 251 251]);

end

%% helper functions
function event_matrix = medi_reshape(vector_in, width, height)
event_matrix = double(reshape(vector_in,width,height));
event_matrix = event_matrix';
end

function corrected_event_matrix = rawhdfcor(event_matrix_in, tot_matrix_in)
event_data_matrix_filter = medfilt2(event_matrix_in, [2 2]); % Define 2x2 Median filter
med_dif = event_matrix_in - event_data_matrix_filter; % compare the original and filtered data sets
med_thresh_pos = 30; % define a threshold difference
med_thresh_neg = -med_thresh_pos;


% apply the median filter
corrected_event_matrix = event_matrix_in;

% For those higher than the upper threshold, or lower than the lower
% threshold, replace the pixel data by the filtered pixel value
corrected_event_matrix(med_dif > med_thresh_pos) = event_data_matrix_filter(med_dif > med_thresh_pos);
corrected_event_matrix(med_dif < med_thresh_neg) = event_data_matrix_filter(med_dif < med_thresh_neg);
% Check the dead pixels based on iTOT data
corrected_event_matrix(tot_matrix_in == 0) = event_data_matrix_filter(tot_matrix_in == 0);
end
