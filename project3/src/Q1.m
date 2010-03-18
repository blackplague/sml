function x = Q1()
%Q1 start
k = 20;
data = load('hands.txt');
%plot (data(k,1:56), data(k, 57:end))

mu = [mean(data(1:40,1)), mean(data(1:40,57))];
sigma = [cov([data(1:40,1),data(1:40,57)])];
%plot_normal_dist(mu, sigma);

mu_vector = zeros(112,1);
hold on
for i=1:56
    
   mu = [mean(data(1:40,i)),mean(data(1:40,i+56))];
   mu_vector(i,1) = mu(1);
   mu_vector(i+56,1) = mu(2);
   sigma = [cov([data(1:40,1),data(1:40,i+56)])];
   %plot_normal_dist(mu,sigma);
end
hold off
size(mu_vector(1:56,1))
size(mu_vector(57:end,1))
%Q1 slut

%Q2 start
%plot mean hand
%plot (mu_vector(1:56,1), mu_vector(57:end,1))

sigma_all = cov(data);
[v,d1] = eig(sigma_all);
d1 = diag(d1);
d = d1/(sum(d1));
e = zeros(size(d)+1,1);
d = flipud(d);
e(2:end,1) = cumsum(d);
%plot(0:112 ,e(1:end,:));
%axis([0 20 0 1]);
hand = 112;
best = v(:, hand)';
%best = best/(sum(best));
d1(hand);
muhaa1 = (mu_vector' - 2*sqrt(d1(hand))*best)';
muhaa2 = (mu_vector' + 2*sqrt(d1(hand))*best)';
hold on
%plot(muhaa1(1:56,1), muhaa1(57:end,1), 'Color', 'green')
%plot(muhaa2(1:56,1), muhaa2(57:end,1), 'Color', 'red')
%plot(mu_vector(1:56,1), mu_vector(57:end,1), 'Color', 'blue')
hold off
%Q2 slut


%Q3 start
k1 = v(:, 112);
k2 = v(:, 111);
lambda1 = d1(112);
lambda2 = d1(111);
us = zeros(40,1);
vs = zeros(40,1);
for i=1:40
    us(i) = 1/(sqrt(lambda1))*k1'*(data(i,:)'-mu_vector);
    vs(i) = 1/(sqrt(lambda2))*k2'*(data(i,:)'-mu_vector);
end
scatter(us,vs)

distences = (us.^2 + vs.^2)
hold on 
%plot (data(40,1:56), data(40, 57:end), 'Color', 'green')
%plot (data(38,1:56), data(38, 57:end), 'Color', 'red')

hold off

%Q3 slut
end