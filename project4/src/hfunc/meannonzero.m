function [ mean_vector ] = meannonzero( path )
%meannonzero returns a mean vector based on non zero entries in each column
%   Detailed explanation goes here

    [~, ~, rv] = getuserinfo( path );
        
    mean_vector = zeros(1,size(rv,2));
    
    for i=1:size(rv,2)
        mean_vector(i) = sum(rv(rv(:,i)>0,i))/length(rv(rv(:,i)>0,1));
    end
end

