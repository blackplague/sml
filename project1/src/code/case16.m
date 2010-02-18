% case 1 Statistical Methods for Machine Learning
% MCMC for multivariate normal 

Sigma = [ 0.3 0.2 ; 0.2 0.2 ] ;
mu = [ 1 ; 1 ];


sampler = 'exact'; % options: {'exact','MH'}
L = 5000; % number of samples
stepsize = 200; % step-size in M-H proposal density
lagmax = 50; % maximum lag in autocorrelation

y = zeros(length(mu),L); % samples
C = chol(Sigma,'lower');

switch sampler
    case 'exact'
        
        y = mu(:,ones(L,1)) + C * randn(length(mu),L) ;
        
    case 'MH'
        
        % calculate inverse for computing the normal density
        % In general, it is more stable to use Cholesky and \ operator
        Sigmainv = inv(Sigma) ; 
        
        % choose initial value
        y(:,1) = - mu ;
        
        for l=2:L
            yprop = y(:,l-1) + stepsize * randn(length(mu),1) ;
            
            % only compute terms that don't cancel in M-H acceptance
            if rand < exp(-0.5*(yprop-mu)'*Sigmainv*(yprop-mu) ...
                          +0.5*(y(:,l-1)-mu)'*Sigmainv*(y(:,l-1)-mu))
            
                y(:,l) = yprop ;
            else
                y(:,l) = y(:,l-1);
            end
        end
       
        
end

% plot contour lines and samples
figure(1)
% con = C * [ cos(0:0.01:2*pi) ; sin(0:0.01:2*pi) ] ; con = con + mu(:,ones(size(con,2),1)) ;
% plot(con(1,:),con(2,:),y(1,:),y(2,:),'*')
dx = -1:0.1:3;
dy = -1:0.1:3;
[ grdx grdy ] = meshgrid( dx, dy );
dens = mvnpdf( [ grdx(:) grdy(:) ], mu', Sigma );
dens = reshape( dens, length( dx ), length( dy ) );
imagesc( dx, dy, dens )
set(gca,'YDir','normal')
hold on
plot(y(1,:),y(2,:),'w*')
hold off

% plot auto-correlation function
figure(2)
muhat = mean(y,2);
sigma2hat = std(y,0,2).^2;

R = zeros(length(mu),lagmax) ;
for k=1:lagmax
    for d=1:length(mu)
        R(d,k) = 1/((L-k)*sigma2hat(d)) * (y(d,1:end-k) - muhat(d) ) *  (y(d,1+k:end) - muhat(d) )' ; 
    end
end

for d=1:length(mu)
    subplot(length(mu),1,d)
    plot(1:k,R(d,:),'o' );
    xlabel( 'lag' ); ylabel( 'corr coeff' )
    grid on
    axis([1 lagmax -1 1])
end
