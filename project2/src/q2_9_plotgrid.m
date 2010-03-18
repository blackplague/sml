clf
x1min=-3;x1max=3;x2min=-3;x2max=3;
[x1,x2]=meshgrid(x1min:0.1:x1max,x2min:0.1:x2max);

p1=conddensity(x1,x2,1);
p2=conddensity(x1,x2,2);

% N1=50;
% N2=250;
% xtr1 = getdataconddensity(N1,1);
% xtr2 = getdataconddensity(N2,2);


figure(1)
contour(x1,x2,p1)
% hold on
% plot(xtr1(:,1),xtr1(:,2),'go','MarkerSize',7)
% hold off

figure(2)
contour(x1,x2,p2)
% hold on
% plot(xtr2(:,1),xtr2(:,2),'ro','MarkerSize',7)
% hold off

figure(3), hold on

pclass1 = zeros(size(x1)); % class probabilities for class 1
pclass2 = zeros(size(x1)); % class probabilities for class 2
S = zeros(size(x1)); % variable needed for edges

loss = [0 1000; 1 0]
% loss = [0 1; 1 0]

[d1,d2] = size(x1);
for i=1:d1
    for j=1:d2
        % insert calculation of class probabilities here!
        pclass1(i,j) = loss(2,1)*p1(i,j) + loss(2,2)*p1(i,j);
        pclass2(i,j) = loss(1,1)*p2(i,j) + loss(1,2)*p2(i,j);
        % insert calculation of class probabilities here!
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