function [ ratedbyno ] = numberofratings( rm )
%numberofratings recieves a movie_id and returns how many who have rated it
%   Detailed explanation goes here

    ratedbyno = cell(1,size(rm,2));

    for i=1:size(rm,2)
        ratedbyno{i} = int2str(length(rm(rm(:,i)>0)));
    end
end

