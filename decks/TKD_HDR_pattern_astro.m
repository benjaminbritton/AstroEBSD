%% Preamble
% TKD_HDR_Pattern.m - a script in AstroEBSD Package to generate a stitched
% high dynamic range on-axis transmission kikuchi diffraction pattern (TKP)
%
% Originally scripted by Tianbi Zhang in November 2022 with major rework in April 2023
% For how this script works, please refer to our paper
% 
% Pre-requisite: MATLAB Image Processing & Statistics and Machine Learning Toolbox
%
% Additional requisites (does not affect this code, but helps you check
% file structure of the h5)
% (a) h5 file viewer to check h5 file structure
%
% Inputs: on-axis transmission Kikuchi patterns captured by a TimePIX3
% direct detector in .h5 format. The patterns can be saved as individual
% frames within a single .h5 file, or as a single integrated frame. There
% are two different functions to import the data to the workspace.
%
% Outputs: stitched high dynamic range pattern (raw and flatfielded), 
% corresponding 2D FFT spectra (modulus) and assorted plots.
%
% The standard version of the script generates similar plots as Figure 3-5
% of the original paper (you will need a simulated pattern)
%
% You can customize this script to produce other plots - please work on a
% copy of this code for your own pattern, and acknowledge the original 
% work if possible.

%%
home;
close all;
clear;

%% Settings
export_option = false; % select true if you want image outputs
multframe_option = false; % true if the data are stored in individual frames within a h5 file
multframe_num = 50; % the number of frames in a multi-frame h5 file. 

% InputUser.Astro_loc='C:\Users\billy\Documents\MATLAB\AstroEBSD_v2-live'; 
% run(fullfile(InputUser.Astro_loc,'start_AstroEBSD.m')); % run AstroEBSD

%% PC info
DD = 10.535; %% CL in mm - you need to measure this in the SEM.
PCx = 136; %%PCx, PCy in pixel units - you can read from the image
PCy = 146; % 149

sat_index = 1022;

%% Import Data

path1= 'C:\Users\tianbi\OneDrive - UBC';
path2 = 'PhD\TKD\Al_HDR\Demo';


% Please fill in the file names in decreasing order of exposure time. The
% file corresponding to the baseline exposure time t0 should be the 
% first one.
h5filenames = {'spot1_50f_01s.h5','spot1_50f_001s.h5', 'spot1_50f_0005s.h5','spot1_50f_0001s_Event.h5'}; 
sampleID =  ["0.1s Exposure", "0.01s Exposure", "0.005s Exposure", "0.001s Exposure"];
[~,numfiles] = size(h5filenames);

rawdata = zeros(254^2,numfiles);
normalizeddata = zeros(254^2,numfiles);
for i=1:numfiles
    filename = fullfile(path1, path2, h5filenames{i});
    if multframe_option == true
        rawdata(:,i) = rawhdfmult(filename, multframe_num);
    else
        if i~=numfiles
        rawdata(:,i) = rawhdfint(filename) ./ multframe_num; % due to data formet issue, use a special function here
        else
           rawdata(:,i) = rawhdfint_special(filename) ./ multframe_num;
        end
    end
end

%% Identify saturated domains 
satdomains = cell(1,numfiles);
for i=1:numfiles
    satdomains{i} = find(rawdata(:,i) >= sat_index);
end

fillingdomains = cell(1,numfiles);
for i=1:numfiles
    if i==1
        fillingdomains{i} = find(rawdata(:,i) < sat_index);
    else
        fillingdomains{i} = setdiff(satdomains{i-1}, satdomains{i});
    end
end

%% Plot raw signal versus scatter angle
[xgrid,ygrid] = meshgrid(1:254,1:254);
angle_matrix = atan(sqrt((xgrid - PCx).^2 + (ygrid - PCy).^2) * 55 / 1000 ./ DD) / pi * 180;

%% Signal Normalization
% Calculate normalization factors
eventratio = zeros(length(fillingdomains{1}),numfiles-1);
for j=2:numfiles
    reference_data = rawdata(:,1);
    current_data = rawdata(:,j);
    eventratio(:,j-1) = reference_data(fillingdomains{1})./ current_data(fillingdomains{1});
end

eventratio(isinf(eventratio)) = NaN;

clear reference_data current_data;

scale_factor = mean(eventratio, 1 ,"omitnan");

for i=1:numfiles
    if i==1
        normalizeddata(:,i) = rawdata(:,i);
    else
        normalizeddata(:,i)= rawdata(:,i) .* scale_factor(i-1);
    end
end

%% Fit reference curve to a smoothing function
notdeadpixel = find(rawdata(:,1) ~=0);
angle_vector = angle_matrix(notdeadpixel);
ref_curve = normalizeddata(:,numfiles);
ref_curve = ref_curve(notdeadpixel);

% set up the fit
[xData, yData] = prepareCurveData( angle_vector, ref_curve );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
excludedPoints = xData < 0.3;
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.8;

% Fit model to data.
[fitresult2, gof] = fit( xData, yData, ft, opts );

%% Stitch the pattern
stitchedpattern = zeros(254^2,1);

for i=1:numfiles
    current_data = normalizeddata(:,i);
    stitchedpattern(fillingdomains{i}) = current_data(fillingdomains{i});
end

stitchedpattern = reshape(stitchedpattern, 254,254);
stitched_mask = medfilt2(stitchedpattern, [3 3]);
stitchedpattern(rawdata(:,1)==0) = stitched_mask(rawdata(:,1)==0);


%% Flat field
normfactor = feval(fitresult2, angle_matrix);
normfactor = reshape(normfactor, [254 254]);
stitchedpattern_flat = stitchedpattern ./ normfactor;

% If your stitched pattern still shows uneven background, consider adding a
% Gaussian filter, for example:
stitched_gaussian = imgaussfilt(stitchedpattern_flat, 15);
stitched_flat_gaus = stitchedpattern_flat ./ stitched_gaussian;

%% Visualization and FFT

stitched_fft = fft2(stitchedpattern);
stitched_flat_fft = fft2(stitchedpattern_flat);
stitched_gaus_fft = fft2(stitched_flat_gaus);

figure;
subplot(2,2,1);
imagesc(stitchedpattern); axis xy; axis image; colormap('gray');axis off;
title("Raw stitched");
subplot(2,2,2);
imagesc(stitched_flat_gaus); axis xy; axis image; colormap('gray');axis off;
title("Flattened stitched");
subplot(2,2,3);
imagesc(log10(abs(fftshift(stitched_fft)))); axis xy; axis image; colormap('gray');axis off;
subplot(2,2,4);
imagesc(log10(abs(fftshift(stitched_gaus_fft)))); axis xy; axis image; colormap('gray');axis off;

%% Assorted Plots

figure;
subplot(1,3,1); %Figure 1: raw data as a function of scattering angle
hold on;
for i=1:numfiles
scatter(angle_matrix(:), rawdata(:,i), 1, 'o', 'filled');
end
grid on;
xlabel("Scatter Angle (degrees)");
ylabel("Event Count");
legend(sampleID);
subplot(1,3,2); %Figure 2: normalized data and the conformity to the reference curve
hold on;
for i=1:numfiles
scatter(angle_matrix(:), normalizeddata(:,i), 1, 'o', 'filled');
end
grid on;
xlabel("Scatter Angle (degrees)");
ylabel("Normalized Event Count");
legend(sampleID);
subplot(1,3,3); %Figure 3: Plot fit with the reference curve
h = plot( fitresult2, xData, yData );
legend( h, 'Reference Curve', 'Smoothing Function Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Scatter angle (degree)', 'Interpreter', 'none' );
ylabel( 'Normalized Event Count', 'Interpreter', 'none' );
grid on;

%% Export images: 16-bit TIFF

stitchedpattern_16 = normalizeto16bit(stitchedpattern);
stitchedpattern_flat_16 = normalizeto16bit(stitchedpattern_flat);
stitched_gaus_16 = normalizeto16bit(stitched_flat_gaus);

% note: below are log10 fft spectrum
stitched_fft_16 = normalizeto16bit(log10(abs(fftshift(stitched_fft))));
stitched_flat_fft_16 = normalizeto16bit(log10(abs(fftshift(stitched_flat_fft))));
stitched_gaus_fft_16 = normalizeto16bit(log10(abs(fftshift(stitched_gaus_fft))));


if export_option == true
    imwrite(uint16(flipud(stitchedpattern_16)), fullfile(path1,path2,'stitched_raw.tif'));
    imwrite(uint16(flipud(stitchedpattern_flat_16)), fullfile(path1,path2,'stitched_flat.tif'));
    imwrite(uint16(flipud(stitched_gaus_16)), fullfile(path1,path2,'stitched_flat_gaus.tif'));
    
    imwrite(uint16(flipud(stitched_fft_16)), fullfile(path1,path2,'stitched_raw_fft.tif'));
    imwrite(uint16(flipud(stitched_flat_fft_16)), fullfile(path1,path2,'stitched_flat_fft.tif'));
    imwrite(uint16(flipud(stitched_gaus_fft_16)), fullfile(path1,path2,'stitched_flat_gaus_fft.tif'));
end

%% Side scripts for other plots
%%% SideScript A: patterns vs exposure time (not normalized)
numpats = numfiles+1;
[size_x, size_y] = size(stitchedpattern);
cropsize = 48;
cropped_rect = [20 20 cropsize-1 cropsize-1];

fft_window = hann(cropsize) * hann(cropsize)';
raw_crop = cell(1, numpats);
raw_crop{end} = imcrop(stitchedpattern, cropped_rect); % raw stitched
sampleID{end+1} = 'Raw stitched';
%this will be exposure time longest to shortest, then stitched

% longsampleID = circshift(sampleID,1);

%% Figure: patterns at different exposure times, short to long, and gaussian
figure;
hold on;
for k=1:numfiles
pattern = reshape(rawdata(:,k),[254 254]);
subplot(1,numpats,k);
imagesc(pattern); axis xy; axis image; axis off; colormap('gray'); colorbar;
title(sampleID{k});
end
subplot(1,numpats, numpats);
imagesc(stitchedpattern); axis xy; axis image; axis off; colormap('gray'); colorbar;
title(sampleID{end});

figure;
hold on;
for k=1:numfiles
pattern = reshape(rawdata(:,k),[254 254]);

gaussian = imgaussfilt(pattern,15);
pattern_gaus = pattern ./ gaussian;
subplot(1,numpats,k);
imagesc(pattern_gaus); axis xy; axis image; axis off; colormap('gray');
end
subplot(1,numpats,numpats);
imagesc(stitched_flat_gaus); axis xy; axis image; axis off; colormap('gray');


%% Side scripts B: extract top rectangle and get fft
for k=1:numfiles
pattern = reshape(rawdata(:,k),[254 254]);
raw_crop{k} = imcrop(pattern, cropped_rect);
end

fft_crop = cell(1,numpats);
fft_crop{end} = log10(abs(fftshift(fft2(raw_crop{end} .* fft_window))));

figure;
for m=1:(numpats)
fft_crop{m} = log10(abs(fftshift(fft2(raw_crop{m} .* fft_window))));
subplot(3,numpats,m);
imagesc(raw_crop{m}); axis xy; axis image; colormap('gray'); axis off;
title(sampleID{m});
subplot(3,numpats,m+numpats);
imagesc(fft_crop{m});  axis xy; axis image; colormap('gray'); axis off;
subplot(3,numpats,m+2*numpats);
imagesc(raw_crop{m} .* fft_window);  axis xy; axis image; colormap('gray'); axis off;
end

% Side script 3: histogram
figure;
histogram(stitchedpattern_16, 200,'EdgeColor','none','Normalization','count');
hold on;
histogram(stitched_gaus_16, 200,'EdgeColor','none','Normalization','count');
xlabel("Normalized grayscale");
ylabel("Pixel count");
legend("Raw stitched","Background corrected");
grid on;

%side script 4:angle contours
f3 = figure;
h2 = axes; 
p2 = imagesc(stitched_flat_gaus);
colormap(h2, 'gray');
axis xy; axis image;axis off;
set(h2,'ydir','normal');
h3 = axes;
[p3,h4] = contour(xgrid,ygrid,angle_matrix,[10 20 30 40 50],'ShowText','on', 'linewidth', 2);
set(h3,'color','none','visible','off','fontsize',22);
axis xy; axis image;
clabel(p3,h4,'color','yellow', 'fontsize',20);
set(h3,'ydir', 'normal');
linkaxes([h2 h3]);

%side script 5: filling zone contours
fillingcontour = zeros(254^2,1);

for i=1:numfiles
    current_data = normalizeddata(:,i);
    fillingcontour(fillingdomains{i}) = i;
end

fillingcontour = reshape(fillingcontour, 254,254);

figure;
contour(fillingcontour);
axis xy; axis image;

%% Side Script 6: line profile across a band - need to plot corrected images
% and relative intensity; almost perfectly perpendicular profiling
newsampleID = sampleID;
newsampleID(end) = "Stitched pattern";
extendedsampleID = newsampleID;
extendedsampleID(end+1) = "Simulated Pattern";

% read the simulated pattern
opts = delimitedTextImportOptions("NumVariables", 254, "Encoding", "UTF-8");
% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "VarName57", "VarName58", "VarName59", "VarName60", "VarName61", "VarName62", "VarName63", "VarName64", "VarName65", "VarName66", "VarName67", "VarName68", "VarName69", "VarName70", "VarName71", "VarName72", "VarName73", "VarName74", "VarName75", "VarName76", "VarName77", "VarName78", "VarName79", "VarName80", "VarName81", "VarName82", "VarName83", "VarName84", "VarName85", "VarName86", "VarName87", "VarName88", "VarName89", "VarName90", "VarName91", "VarName92", "VarName93", "VarName94", "VarName95", "VarName96", "VarName97", "VarName98", "VarName99", "VarName100", "VarName101", "VarName102", "VarName103", "VarName104", "VarName105", "VarName106", "VarName107", "VarName108", "VarName109", "VarName110", "VarName111", "VarName112", "VarName113", "VarName114", "VarName115", "VarName116", "VarName117", "VarName118", "VarName119", "VarName120", "VarName121", "VarName122", "VarName123", "VarName124", "VarName125", "VarName126", "VarName127", "VarName128", "VarName129", "VarName130", "VarName131", "VarName132", "VarName133", "VarName134", "VarName135", "VarName136", "VarName137", "VarName138", "VarName139", "VarName140", "VarName141", "VarName142", "VarName143", "VarName144", "VarName145", "VarName146", "VarName147", "VarName148", "VarName149", "VarName150", "VarName151", "VarName152", "VarName153", "VarName154", "VarName155", "VarName156", "VarName157", "VarName158", "VarName159", "VarName160", "VarName161", "VarName162", "VarName163", "VarName164", "VarName165", "VarName166", "VarName167", "VarName168", "VarName169", "VarName170", "VarName171", "VarName172", "VarName173", "VarName174", "VarName175", "VarName176", "VarName177", "VarName178", "VarName179", "VarName180", "VarName181", "VarName182", "VarName183", "VarName184", "VarName185", "VarName186", "VarName187", "VarName188", "VarName189", "VarName190", "VarName191", "VarName192", "VarName193", "VarName194", "VarName195", "VarName196", "VarName197", "VarName198", "VarName199", "VarName200", "VarName201", "VarName202", "VarName203", "VarName204", "VarName205", "VarName206", "VarName207", "VarName208", "VarName209", "VarName210", "VarName211", "VarName212", "VarName213", "VarName214", "VarName215", "VarName216", "VarName217", "VarName218", "VarName219", "VarName220", "VarName221", "VarName222", "VarName223", "VarName224", "VarName225", "VarName226", "VarName227", "VarName228", "VarName229", "VarName230", "VarName231", "VarName232", "VarName233", "VarName234", "VarName235", "VarName236", "VarName237", "VarName238", "VarName239", "VarName240", "VarName241", "VarName242", "VarName243", "VarName244", "VarName245", "VarName246", "VarName247", "VarName248", "VarName249", "VarName250", "VarName251", "VarName252", "VarName253", "VarName254"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
simpattern_table = fullfile(path1, path2, 'sim_pattern.csv');

simpattern = readtable(simpattern_table, opts);
simpattern = table2array(simpattern);

% meansim = mean(simpattern(:));
% meanstitched = mean(stitched_flat_gaus(:));

simpattern = simpattern + 0.5;
clear opts;

% simpattern = imresize(simpattern, [254 254]);
figure;
imagesc(simpattern); axis xy; axis image; colormap('gray'); axis off;

[line_1_x, line_1_y] = fitedge(28, 57, 34, 31); % these are read from the pattern manually
[line_2_x, line_2_y] = fitedge(200,222,167,178);

figure;
subplot(2,2,[1 3]);
hold on;
imagesc(stitched_flat_gaus); axis xy; axis image; colormap('gray'); axis off;
hold on;
plot(line_1_x, line_1_y,'Color','cyan','LineWidth',2); 
plot(line_2_x, line_2_y,'Color','magenta','LineWidth',2); 
hold off;


for k=1:numfiles
thispattern = reshape(normalizeddata(:,k),[254 254]);
thisnormfactor = feval(fitresult2, angle_matrix(:));
thisnormfactor = reshape(thisnormfactor,[254 254]);
bgcorpattern = thispattern ./ thisnormfactor;
subplot(2,2,2);
hold on;
lineprof1 = improfile(bgcorpattern, line_1_x, line_1_y, 'nearest');
plot(lineprof1,'LineWidth',2);
subplot(2,2,4);
hold on;
lineprof2 = improfile(bgcorpattern, line_2_x, line_2_y, 'nearest');
plot(lineprof2,'LineWidth',2);
end
subplot(2,2,2);
lineprof = improfile(stitchedpattern_flat, line_1_x,line_1_y, 'nearest');
plot(lineprof,'LineWidth',2);
lineprof = improfile(simpattern, line_1_x,line_1_y, 'nearest');
plot(lineprof,'LineWidth',2);
legend(extendedsampleID);
grid on;
xlabel("Pixel Position");
ylabel("Corrected Intensity (a.u.)");
title("Line 1 (Cyan)");
subplot(2,2,4);
lineprof = improfile(stitchedpattern_flat, line_2_x, line_2_y, 'nearest');
plot(lineprof,'LineWidth',2);
lineprof = improfile(simpattern, line_2_x-2,line_2_y, 'nearest');
plot(lineprof,'LineWidth',2);
legend(extendedsampleID);
grid on;
xlabel("Pixel Position");
ylabel("Corrected Intensity (a.u.)");
title("Line 2 (Magenta)");

%% Helper functions
function event_data_vector = rawhdfint(h5filename)

% Read width and height values from the HDF5 file. Both should be 256
% because of the DED geometry
width = h5read(h5filename,'/Frame_0/Width');
height = h5read(h5filename,'/Frame_0/Height');

% Extract pixel-by-pixel intensity from the "Event" and "iTOT" data.
% Change "Event" to "iToT", "ToT", "ToA" etc. based on data structure.
% Please check the file under HDF5 viewer or use h5info() first.
event_data = h5read(h5filename,'/Frame_0//SubFrames/Event/Data');
event_data_matrix = double(reshape(event_data,width,height));

event_data_matrix = event_data_matrix';
event_data_matrix = imcrop(event_data_matrix, [2 2 253 253]);

event_data_matrix(11,60) = 0;
event_data_matrix(26,208) = 0;
event_data_matrix(240,22) = 0;
event_data_matrix(44,47) = 0;

event_data_vector = double(reshape(event_data_matrix, [],1));
end

function event_data_vector = rawhdfint_special(h5filename)
% Read width and height values from the HDF5 file. Both should be 256
% because of the DED geometry
width = h5read(h5filename,'/Frame_0/Width');
height = h5read(h5filename,'/Frame_0/Height');

% Extract pixel-by-pixel intensity from the "Event" and "iTOT" data.
% Change "Event" to "iToT", "ToT", "ToA" etc. based on data structure.
% Please check the file under HDF5 viewer or use h5info() first.
event_data = h5read(h5filename,'/Frame_0/Data');
event_data_matrix = double(reshape(event_data,width,height));

event_data_matrix = event_data_matrix';
event_data_matrix = imcrop(event_data_matrix, [2 2 253 253]);

event_data_matrix(11,60) = 0;
event_data_matrix(26,208) = 0;
event_data_matrix(240,22) = 0;
event_data_matrix(44,47) = 0;

event_data_vector = double(reshape(event_data_matrix, [],1));
end


function event_data_vector = rawhdfmult(h5filename, framenum)
width = h5read(h5filename,'/Frame_0/Width');
height = h5read(h5filename,'/Frame_0/Height');

event_data = double(zeros(65536,1));

for i=0:(framenum-1)
    frame_event_key = strcat('/Frame_',num2str(i),'/SubFrames/Event/Data');
    event_data_this_frame = double(h5read(h5filename, frame_event_key));
    
    event_data = event_data + event_data_this_frame;
end

event_data_matrix = double(reshape(event_data,width,height));

event_data_matrix = event_data_matrix';
event_data_matrix = imcrop(event_data_matrix, [2 2 253 253]);

event_data_matrix(11,60) = 0;
event_data_matrix(26,208) = 0;
event_data_matrix(240,22) = 0;
event_data_matrix(44,47) = 0;

event_data_matrix = event_data_matrix / framenum;

event_data_vector = double(reshape(event_data_matrix, [],1));
end

% write a dead pixel correction program

function normalized_matrix = normalizeto1(matrix_in)
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
end

function normalized_matrix = normalizeto16bit(matrix_in)
% normalized_matrix = matrix_in;
normalized_matrix = matrix_in - min(matrix_in(:)); % make the lowest 0
normalized_matrix = normalized_matrix ./ max(normalized_matrix(:)); % Normalize
normalized_matrix = uint16(normalized_matrix * (2^16 - 1)); % Convert to 16 bit
end

function [fitlinex, fitliney] = fitedge(x1,x2,y1,y2)
fitslope = (y2 - y1) / (x2 - x1);

fitlinex = linspace(x1,x2,(abs(x2 - x1) + 1));
fitliney = y1 + fitslope .* (fitlinex - x1);
end
