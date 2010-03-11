function [ xout n ] = Normlized_Histogram( data, bins )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n xout] = hist(data, bins);

bin_width = abs(xout(1) - xout(2)); % Bin width

area = sum(n)*bin_width;

n = n / area;

bar(xout, n, 'BarWidth', 1)

end

