%% hdfbatchebsd.m
% This script performs background subtraction afor a batch of input hdf
% patterns one by one and export them to images (.tif by default).
% Assumes EBSD grometry!

%% Toolbox locations for AstroEBSD and MTEX
InputUser.Astro_loc='C:\Users\tianbi\Documents\MATLAB\AstroEBSD'; 

% Initialize AstroEBSD and MTEX
% run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m'));

%% HDF5 file I/O
path = 'C:\Users\billy\OneDrive - UBC\PhD\TKD\20230807_Cu_EBSD';
file_pre = 'Spot';

maxindex = 400;
opts.output = 1;
gauss = 10; % sigma for Gaussian filter
gamma = 1; % Gamma filter (1 by default which means no change)

for i=1:maxindex
    i_Padded = int2str(i);
    filename = strcat(file_pre, i_Padded);

    [raw_pattern, bgcor_pattern] = hdfpatternprocess(path, filename, 25, 0.9, opts.output, 1);
    
    fprintf('%d of %d done.\n',i,maxindex);

end
