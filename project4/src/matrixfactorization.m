function [Y, U, V, RMSE] = matrixfactorization(X,K,lambda, eta, tmax)
% [Y,U,V] = matrixfactorization(X,K,lambda)
%
%  Y = U * V' - predictions on all entries in X
%  U - user factor matrix
%  V - item matrix
%  X - User-Item (Movie) matrix
%  K - the dimensionality of the latent space (default 10)
%  lambda - weight decay parameter (default 0)

if nargin < 2
    K = 10;
end

if nargin < 3
    lambda = 0 ;
end

if nargin < 4
    eta = 0.0005; % learning rate - rather conservtive more seems necessary to get convergence
end

if nargin < 5
    tmax = 500; % maximum number of iterations
end
    
[N,M] = size(X);

indx = X > 0 ; % find ranked movies - will be used for logical indexing
Ntr = sum(sum(indx)); % number of traning examples

U = 0.1*randn(N,K); % initialize to small random values
V = 0.1*randn(M,K); % initialize to small random values

Error_c = inf;
for t=1:tmax

    E = X - U*V' ; % compute residual for all entries
    
    RMSE = sqrt(sum(sum((indx .* E).^2))/Ntr); % compute fit on training data 
    
%     % make step size adaptive
%     if Error<Error_c
%         eta = 1.05*eta;
%     else
%         eta = 0.5*eta;
%     end
%     
    % loop over users
    for n=1:N
        U(n,:) =  U(n,:) + eta * (E(n,indx(n,:)) * V(indx(n,:),:) - lambda * U(n,:));  
    end
        
    E = X - U*V' ; % compute residual for all entries
    
    % loop over items
    for m=1:M
        V(m,:) =  V(m,:) + eta * (E(indx(:,m),m)' * U(indx(:,m),:) - lambda * V(m,:));   
    end
end
Y = U*V';
