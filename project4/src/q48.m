function G = q48( path )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    function G = googleMatrix( path )
    %googleMatrix - Constructs the google matrix having the size
    %users+movies x users+movies

        % Read in the rating matrix R given the path.
        [ R, ~, ~, ~, ~ ] = ...
            readmovielens('ml-data/', path );

        % Init sizes.
        userSize = size(R,1); movieSize = size(R,2);
        % Initialize G to zeros.
        G = zeros(userSize+movieSize,userSize+movieSize);
        
        % Fill in ones where the R matrix entries are r > 0
        G(1:userSize,userSize+1:end) = R>0;
        G(userSize+1:end,1:userSize) = (R>0)';
    end

if nargin < 1
    path = 'u.data';
end

    function G = run( path )
        G = googleMatrix( path );
    end
    G = run( path );
end