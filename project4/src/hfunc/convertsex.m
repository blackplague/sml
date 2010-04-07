function [sexmatrix] = convertsex( sexes )
%convertsex Function to convert the gender strings to numerical values,
% it maps F -> 1 and M -> 0.

    sexmatrix = zeros(size(sexes));

    for i = 1:size(sexes,1)
        if strcmp('F', sexes{i})
            sexmatrix(i,1) = 1;
        end
    end
end