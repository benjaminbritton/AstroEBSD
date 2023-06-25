function [rawdata,numfiles] = h5_pixet(h5filenames,path1,multframe_option,multframe_num)
%H5_PIXET read patterns from a h5 pixet output
% (c) Tianbi Zhang and Ben Britton, 2023

[~,numfiles] = size(h5filenames);

rawdata = zeros(254^2,numfiles);

for i=1:numfiles
    filename = fullfile(path1, h5filenames{i});
    if exist(filename,'file') == false
        error(['The input file ' filename ' cannot be found']);
    end

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

end

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

%TZ QUERY
%what are these pixels?
pix_cleanxval=[11,26,240,44];
pix_cleanyval=[60,208,22,47];
event_data_matrix=pix_clean(event_data_matrix,pix_cleanxval,pix_cleanyval);

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

pix_cleanxval=[11,26,240,44];
pix_cleanyval=[60,208,22,47];
event_data_matrix=pix_clean(event_data_matrix,pix_cleanxval,pix_cleanyval);

event_data_vector = double(reshape(event_data_matrix, [],1));
end

function h5matrix=pix_clean(h5matrix,xval,yval)
%clean specific pixels fro the array to 0
for n=1:numel(xval)
    h5matrix(xval(n),yval(n)) = 0;
end

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

%TZ QUERY
pix_cleanxval=[11,26,240,44];
pix_cleanyval=[60,208,22,47];
event_data_matrix=pix_clean(event_data_matrix,pix_cleanxval,pix_cleanyval);


event_data_matrix = event_data_matrix / framenum;

event_data_vector = double(reshape(event_data_matrix, [],1));
end