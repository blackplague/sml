function T = createTransitionMatrix()
%createTransitionMatrix - returns T_ij = G_ij/L_j, where
% L_j = Sum_i(G_ij)

    % Get Google matrix from q48.
    G = q48();
        
    T = zeros(size(G));
    L = sum(G);
        
    for i=1:size(G,1)
        T(:,i) = G(:,i)/L(i);
    end
end