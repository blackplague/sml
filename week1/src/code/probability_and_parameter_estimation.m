%Assinment 1
function result = probability_and_parameter_estimation()
sigma = [0.3 0.2 ; 0.2 0.2];
mu = [1.0 1.0];
N = 100;
x = Q1(sigma, mu, N);
%Q2(mu, x, N);
Q3(x);
%Q4(x);
Q5();
result = 'indet lige nu'
end

%Questian 1.1

function x = Q1(sigma, mu, N)
L = chol(sigma);
z = randn(N,2);
x = repmat(mu,N,1) + z*L;
scatter(x(:,1),x(:,2));
end

%questian 1.2
function muml = Q2(mu,x, N)
muml = sum(x)/N;
con = sum(x-repmat(muml,N,1))/(N-1);
quantify = sqrt((mu(1)-muml(1))^2+(mu(2)-muml(2))^2);

%scatter(muml(1),numl(2));
end

%questian 1.3
function k = Q3(x)
nbin = 13;
hist(x(:,1),nbin);
%hist(x(:,2),nbin);
k = 0;
end

%questian 1.4
function k = Q4(x)
nbin = [20, 20];
hist3(x, nbin);
k = 0;
end

%questian 1.5
function k = Q5()
absdif = zeros(3,1);
lambda = 1.5;
mu_y = 1/lambda;
for p = 1:3
    L = 10  
    L = L^p;
    y = -lambda^(-1)*log(1-rand(L,1000));
    yhat = mean(mean(y,1),2);
    absdif(p,1) = abs(mu_y - yhat);
end
plot([10;100;1000],absdif)
k = 0;
end