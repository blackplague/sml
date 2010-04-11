function nearest = findnearest(nonearest, movie_id, U1, U2)

    distance = zeros(length(U1),length(U2));

%     distance = pdist(U);
    
    % Find the distance between every point and generate a symmetric
    % distance matrix.
    for i = 1:length(U1)
        for j = 1:length(U2)
            distance(i,j) = sqrt((U1(i)-U1(j))^2+(U2(i)-U2(j))^2);
        end
    end

    findin = distance(movie_id,:);

    [~, IDX] = sort(findin);

    nearest = IDX(1:nonearest+1);
    
% Own original code, ejects the one we search out from, and return the ten
% next ones.
%     nearest = zeros(nonearest,1);
%     % Eject movie_id
%     [~, I] = min(findin);
%     findin(I) = max(findin)+10e6;
%     
%     for i = 1:length(nearest)
%         [~, I] = min(findin);
%         nearest(i) = I;
%         findin(I) = max(findin)+10e6;
%     end
end