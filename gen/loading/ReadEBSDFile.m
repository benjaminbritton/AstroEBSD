function [pattern1] = ReadEBSDFile( patternfile,fliptype )
%PatternReadFile Reads an image file from disk
%[pattern1] = ReadEBSDFile( patternfile,fliptype )
%
%INPUTS 
% patternfile = file path
% fliptype = 1 = UD; 2 = LR; 3 = UD + LR (optional)
%
%OUTPUTS
% pattern1 = output pattern

if exist('fliptype','var') == 0
    fliptype = 0;
end

if exist(patternfile,'file') ~= 2
    error(['The pattern file ' patternfile ' does not exist']);
end

%read the pattern
pattern1=imread(patternfile);

%check if it is a RGB file
if size(pattern1,3) == 3
    pattern1=rgb2gray(pattern1);
end

%convert to double
if isa(pattern1,'double') ~=1
pattern1=im2double(pattern1);
end

if fliptype == 1 || fliptype == 3%up
    pattern1=flipud(pattern1);
end
    if fliptype == 2 || fliptype == 3%lr
        pattern1=fliplr(pattern1);
    end
    
end

