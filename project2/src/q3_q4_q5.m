clf
clear
close all
clc
x1min=-3;x1max=3;x2min=-3;x2max=3;
[x1,x2]= meshgrid(x1min:0.1:x1max,x2min:0.1:x2max);

p1=conddensity(x1,x2,1);
p2=conddensity(x1,x2,2);

N1=100;
N2=100;

xtr1 = getdataconddensity(N1,1);
xtr2 = getdataconddensity(N2,2);

% start på opve 2.4
mu1 = mean(xtr1', 2);
mu2 = mean(xtr2', 2);

S1 = cov(xtr1);
S2 = cov(xtr2);

S = (N1/(N1+N2))*S1+(N2/(N1+N2))*S2;

destibution1 = zeros(size(p1));
[d1,d2] = size(p1);
Pi = (N1/(N1+N2));
 for i=1:d1
     for j=1:d2
     destribution1(i,j) = Pi*1/(2*pi*sqrt(det(S)))*exp(-0.5* ...
             ([x1(i,j), x2(i,j)]-[mu1(1),mu1(2)])*...
             inv(S)*([x1(i,j), x2(i,j)]-[mu1(1),mu1(2)])');
     end
 end
 
 

destibution2 = zeros(size(p2));
[d1,d2] = size(p2);

 for i=1:d1
     for j=1:d2
     destribution2(i,j) = (1-Pi)*1/(2*pi*sqrt(det(S)))*exp(-0.5* ...
             ([x1(i,j), x2(i,j)]-[mu2(1),mu2(2)])*...
             inv(S)*([x1(i,j), x2(i,j)]-[mu2(1),mu2(2)])');
     end
 end
% slut på opgave 2.4
 
% start på opgave 2.5
t = [zeros(1,length(xtr1)),ones(1,length(xtr2))];

X = [xtr1; xtr2];
size(X);

Phi1 = [ones(length(X),1),X];
Phi2 = [ones(length(X),1),X,X.*X];
size(Phi2)
% w_old = [0.0001,0.0001,0.0001, 0.0001, 0.0001]'; % til Phi2
w_old = [0.0001,0.0001,0.0001]';
a_n = Phi1 * w_old;
y = zeros(size(a_n));

for i=1:length(y)
   y(i) = 1/(1+exp(-a_n(i))); 
end


R = zeros(length(y),length(y));

for i=1:length(y)
    R(i,i) = y(i)*(1-y(i));
end
epcilon = 400;
u = 14;
% til Phi2 endre alle til Phi2
while(u > 1)
    w_new = inv(Phi1'*R*Phi1)*Phi1'*R*(Phi1*w_old-inv(R)*(y-t'));
    sum(sum(w_new))-sum(sum(w_old));
    if (sum(sum(w_new))-sum(sum(w_old))< epcilon)
        u = u - 1;
    end
    w_old = w_new;
end
a_n = Phi1*w_new;

% slut på opgave 2.5



figure(1)
contour(x1,x2,p1)
%hold on
%plot(xtr1(:,1),xtr1(:,2),'go','MarkerSize',7)
%hold off

figure(2)
contour(x1,x2,p2)
% hold on
% plot(xtr2(:,1),xtr2(:,2),'ro','MarkerSize',7)
% hold off

figure(3), hold on
pclass1 = zeros(size(x1)); % class probabilities for class 1
pclass2 = zeros(size(x1)); % class probabilities for class 2
S = zeros(size(x1)); % variable needed for edges
%opgave 2.5
P = zeros(size(x1));
[d1,d2] = size(x1);
for i=1:d1
    for j=1:d2
        x = [x1(i,j), x2(i,j)];
        phi = [1 x]';
        if((1/(1+exp(-(w_new'*phi))))>0.5)
            pclass1(i,j) = 1;
            pclass1(i,j) = 0;
        else
            pclass1(i,j) = 0;
            pclass1(i,j) = 1;
        end
    end
end
%opgave 2.5 slut

sump1 = sum(sum(p1));
sump2 = sum(sum(p2));
sumc1 = sump1/sump2;
sumc2 = sump2/sump1;

[d1,d2] = size(x1);
for i=1:d1
    for j=1:d2
        %opgave 2.5 ingen markeret
        % Orginal
        %pclass1(i,j)=p1(i,j);
        %pclass2(i,j)=p2(i,j);
        % Bayes 2.3
        %pclass1(i,j)=p1(i,j)*sumc1/(1/(60*60)); 
        %pclass2(i,j)=p2(i,j)*sumc2/(1/(60*60));
        % Maximum likelihood 2.4
        %pclass1(i,j)=destribution1(i,j);
        %pclass2(i,j)=destribution2(i,j);
        if  pclass1(i,j) > pclass2(i,j)  % plot grid point as green
            plot(x1(i,j),x2(i,j),'g*','MarkerSize',2)
            S(i,j)=1;
        else % plot grid point as red
            plot(x1(i,j),x2(i,j),'r*','MarkerSize',2)
            S(i,j)=0;
        end
    end
end

% plot(xtr1(:,1),xtr1(:,2),'go','MarkerSize',7)
% plot(xtr2(:,1),xtr2(:,2),'ro','MarkerSize',7)

axis([x1min x1max x2min x2max])
hold off

% % code to generate edges
% E = edge(S); % find edge
% [Ex, Ey] = find(E); % find x and y indices of edge
% %plot(x1(Ex),x2(Ey),'.', 'MarkerSize',4.5) - plot the edge detection
%  
% xx1 = x1(1,:)';
% xx2 = x2(:,1);
% [xe(1),Ix] = min(xx1(Ex)); % starting point for edge
% Ex(Ix)=[];
% ye(1)=xx2(Ey(Ix)); 
% Ey(Ix)=[];
% h=abs(xx1(1)-xx1(2)); % edge length between two grid points 
% i=2;
% while ~isempty(Ex)
%     [dist,Id]=sort(sqrt((xx2(Ey)-ye(i-1)).^2+(xx1(Ex)-xe(i-1)).^2)); % calculate nearest point
%     if dist(1)<2*h % if closer than 2 times the edge length: connect to previous point
%         xe(i)=xx1(Ex(Id(1))); 
%         Ex(Id(1))=[];
%         ye(i)=x2(Ey(Id(1)));
%         Ey(Id(1))=[];
%     else % if farther away: make new line and plot old line
%         y_smooth=(ye(1:end-2)+ye(2:end-1)+ye(3:end))/3; % smooth y - this could be optimized
%         plot(xe, [ye(1), y_smooth, ye(end)],'-k') % plot edge as 'smooth' function 
%         xe=[]; ye=[]; i=1;
%         [xe(i),Ix] = min(xx1(Ex)); % starting point
%         Ex(Ix)=[];
%         ye(i)=xx2(Ey(Ix)); 
%         Ey(Ix)=[];
%     end
%     i=i+1;
% end
% y_smooth=(ye(1:end-2)+ye(2:end-1)+ye(3:end))/3; % smooth y - this could be optimized
% %plot(xe, [ye(1), y_smooth, ye(end)],'-k') % plot edge as 'smooth' function 
% plot([ye(1), y_smooth, ye(end)],xe,'-k') % plot edge as 'smooth' function 
% axis([x1min x1max x2min x2max])
% hold off
