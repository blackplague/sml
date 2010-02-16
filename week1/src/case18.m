%% Read Image
RGB = im2double (imread ('./images/kande1.jpg'));

%% Define training and test regions
tll = [328, 150]; % training: lower left
tur = [264, 330]; % training: upper right

%% Show image and regions
imshow (RGB);
hold on
  plot ([tur(2), tll(2), tll(2), tur(2), tur(2)], ...
        [tur(1), tur(1), tll(1), tll(1), tur(1)], ...
        'k-','linewidth', 2,'Color', 'green');
  text (tll(2), tur(1)-15, 'Training');

hold off

%% Extract regions
training_data = RGB (tur(1):tll(1), tll(2):tur(2), :);
training_data = reshape (training_data, size (training_data, 1) * size (training_data, 2), 3);

RGB = im2double (imread ('./images/kande1.jpg'));
%RGB now in only one vector
reshape_data = reshape (RGB, size (RGB, 1) * size (RGB, 2), 3);

%calculate the mean an coveriance of the trainings data
means = mean(training_data)
coveriancen = inv(cov(training_data))

%use the mean an coveriance in are gauss detribution on RGB
Z = size(reshape_data, 1)
%P = zeros(Z, 3);

tmp = 1/((2*pi)^(3/2)*det(coveriancen)^(1/2));

xmu = reshape_data-repmat(means,Z,1);

%calculate the coveriancen^-1*(x-mu) in vector form
sigmaxmu1 = coveriancen(1,1)*xmu(:,1)+coveriancen(1,2)*xmu(:,2)+coveriancen(1,3)*xmu(:,3);
sigmaxmu2 = coveriancen(2,1)*xmu(:,1)+coveriancen(2,2)*xmu(:,2)+coveriancen(2,3)*xmu(:,3);
sigmaxmu3 = coveriancen(3,1)*xmu(:,1)+coveriancen(3,2)*xmu(:,2)+coveriancen(3,3)*xmu(:,3);

sigma_tims_data = [sigmaxmu1,sigmaxmu2,sigmaxmu3];
size(sigma_tims_data)
p = zeros(Z,1);

for i = 1:Z
   p(i,1) = tmp*exp(-0.5*(xmu(i,:))*sigma_tims_data(i,:)');
end

summ = sum(p);
new_image = reshape (p,  size (RGB, 1), size (RGB, 2),1);
%imagesc(-new_image), colormap(gray);


pos = zeros(Z,2);
q_hat = [0,0];


for i = 1:Z
    q_hat = q_hat+p(i,:)*([floor((i-1)/(size(RGB,1)))+1,mod((i-1),size(RGB,1))]);
end 

q_hat = round(q_hat/sum(p));
RGB(q_hat(2)-2:q_hat(2)+2,q_hat(1)-2:q_hat(1)+2,:) = 0;
%imshow (RGB);


([1,2]-q_hat)'
([1,2]-q_hat)

C = [0,0;0,0];
for i = 1:size(new_image,1)
    for j = 1:size(new_image,2)
        C = C + (([i,j]-q_hat)'*([i,j]-q_hat))*new_image(i,j);
    end
end
C = C/(sum(p))
plot_results(RGB,q_hat',C)

%hold on
%  plot(q_hat(1),q_hat(2));%,'k-','linewidth', 2,'Color', 'green');
%hold off

%hold on
%    plot ([tur(2), tll(2), tll(2), tur(2), tur(2)], ...
%        [tur(1), tur(1), tll(1), tll(1), tur(1)], ...
%        'Color','red','k-', 'linewidth', 2);
%hold off
%q_hat = 1/(pixels)*sum(p*)

%lej med at lave den unden loops
%size(tmp)
%size(xmu)
%xmu.*tmp
%P = 1/((2*pi^(3/2)*det(coveriancen)^(1/2)))*
%P = (xmu).*tmp;     dette vil give en matrice hvor de digunale ellemeter i
%matricen er de værdier vi vil have
%size(p)
%new_image = reshape (p,  size (RGB, 1), size (RGB, 2),3);
