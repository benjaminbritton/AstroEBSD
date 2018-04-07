function [ text_out ] = pNameClock
%PNAMECLOCK Summary of this function goes here
%   Detailed explanation goes here

t2=clock;

text_out=[sprintf('%04.0f',t2(1)) '_' sprintf('%02.0f',t2(2)) '_' sprintf('%02.0f',t2(3))...
    '_' sprintf('%02.0f',t2(4)) '_' sprintf('%02.0f',t2(5)) '_' sprintf('%02.0f',t2(6))];


end

