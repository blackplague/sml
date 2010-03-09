function [a] = assinment_2()
%questian 1,
a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;

for i=1:50
    [train, test] = readbodyfat;
    %I Q1 bliver RMS1 lavest så den er best
    [Phi1, Phi2, Phi1_test, Phi2_test, train_data, test_data, a1, b1, c1, d1] = Q1(train, test);
    a = [a,a1];
    b = [b,b1];
    c = [c,c1];
    d = [d,d1];
end
%the result for the Q1
RMS1_train = sum(a)/50;
RMS1_train_max = max(a(2:end));
RMS1_train_min = min(a(2:end));
RMS1_test = sum(b)/50;
RMS1_test_max = max(b(2:end));
RMS1_test_min =min(b(2:end));
RMS2_trainsum = sum(c)/50;
RMS2_train_max = max(c(2:end));
RMS2_train_min = min(c(2:end));
RMS2_test = sum(d)/50;
RMS2_test_max = max(d(2:end));
RMS2_test_min = min(d(2:end));
% questian 2

Q2(Phi1, Phi2, Phi1_test, Phi2_test);
a = 0;
end

function [Phi1, Phi2, Phi1_test, Phi2_test, train_data, test_data, RMS1_train, RMS1_test, RMS2_train, RMS2_test] = Q1(train, test)
%initilice phi1
Phi1 = [ones(1,200)', train(:,4), train(:,7), train(:,8), train(:,9)];
Phi1_test = [ones(1,52)', test(:,4), test(:,7), test(:,8), test(:,9)];

train_data = train(:,2);
test_data = test(:,2);

w_ml_1 = inv(Phi1'*Phi1)*Phi1'*train_data;
y_1 = Phi1 * w_ml_1;
y_1_test = Phi1_test * w_ml_1;

RMS1_train = sqrt(sum((train_data-y_1).^2)/200);
RMS1_test =  sqrt(sum((test_data-y_1_test).^2)/52);

Phi2 = [ones(1,200)', train(:,8)];
Phi2_test = [ones(1,52)', test(:,8)];

w_ml_2 = inv(Phi2'*Phi2)*Phi2'*train_data;
y_2 = Phi2 *w_ml_2;
y_2_test = Phi2_test*w_ml_2;

RMS2_train = sqrt(sum((train_data-y_2).^2)/200);
RMS2_test = sqrt(sum((test_data-y_2_test).^2)/52);


end

function [a] = Q2(Phi1, Phi2, Phi1_test, Phi2_test)
hold on
for j=1:1
    j
[train, test] = readbodyfat;
train_data = train(:,2);
test_data = test(:,2);
beta = 1;
alpha = 1;

for i=1:200
alpha = i;
inv_Sn_1 = alpha*eye(5)+beta*Phi1'*Phi1;
%Detter er vores p(w|t) mean
Mn_1 = beta*inv(inv_Sn_1)*Phi1'*train_data;
y_1 = mean(abs(Phi1*Mn_1-train_data));


RMS1 = sqrt(mean((train_data - Phi1*Mn_1).^2));


inv_Sn_1_test = alpha*eye(5)+beta*Phi1_test'*Phi1_test;
Mn_1_test = beta*inv(inv_Sn_1_test)*Phi1_test'*test_data;
y_1_test = mean(abs(Phi1_test*Mn_1_test-test_data));

RMS2 = sqrt(mean((test_data - Phi1_test*Mn_1_test).^2));


inv_Sn_2 = alpha*eye(1)+beta*Phi2'*Phi2;
Mn_2 = beta*inv(inv_Sn_2)*Phi2'*train_data;
y_2 = mean(abs(Phi2*Mn_2));

RMS3 = sqrt(mean((train_data - Phi2*Mn_2).^2));

inv_Sn_2_test = alpha*eye(1)+beta*Phi2_test'*Phi2_test;
Mn_2_test = beta*inv(inv_Sn_2_test)*Phi2_test'*test_data;
y_2_test = mean(abs(Phi2_test*Mn_2_test-test_data));

RMS4 = sqrt(mean((test_data - Phi2_test*Mn_2_test).^2));

plot(i,RMS1, 'Color', 'green');
plot(i,RMS2, 'Color', 'yellow');
plot(i,RMS3, 'Color', 'red');
plot(i,RMS4, 'Color', 'blue');
end
legend( 'RMS1Traning', 'RMS1Test', 'RMS2Traning', 'RMS2Test', 'Location', 'NorthWest' )

end
hold off

a = 0;

end