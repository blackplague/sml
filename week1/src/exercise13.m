function [ output_args ] = exercise13( y, bins )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x1 = y(:,1);
    x2 = y(:,2);

    subplot(1,2,1)
    hist(x1, bins)
    subplot(1,2,2)
    hist(x2, bins)
    
end

