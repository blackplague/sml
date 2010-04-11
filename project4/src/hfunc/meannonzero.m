function [ mean_vector ] = meannonzero( rm )
%meannonzero returns a mean vector based on non zero entries in each column
%   Detailed explanation goes here
        
    mean_vector = zeros(1,size(rm,2));
    
    for i=1:size(rm,2)
        if (~isempty(rm(rm(:,i)>0,1)))
            mean_vector(i) = sum(rm(rm(:,i)>0,i))/length(rm(rm(:,i)>0,1));
        else
            mean_vector(i) = 0;
        end
    end
end

