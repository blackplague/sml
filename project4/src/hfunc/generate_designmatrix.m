function designmatrix = generate_designmatrix(selection,phifn)
%generate_designmatrix used to generate a design matrix given the function
%handle phifn, e.g. phifn = @(x)[1 x x.^2]
    sphi = size(phifn(1),2)-1;
    M = size(selection,1);
    N = (size(selection,2)*sphi)+1;
    designmatrix = zeros(M,N);
    for row=1:1:M
        designmatrix(row,:) = phifn(selection(row,:))';
    end
end