function [ sigma_est ] = exercise12b( mu_ml, y )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    y = y';

    N = length(y);
    
    temp = ones(100,2);
    
    tmp = [mu_ml(1) * temp(:,1)'; mu_ml(2) * temp(:,2)'];
   
    sigma_est = zeros(2,2);
    
    sigma_est = 1 / (N - 1) * (y - tmp)*(y - tmp)';

end

