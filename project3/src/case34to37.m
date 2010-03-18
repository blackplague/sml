% reading in the tab separeted file with header
rawdata = importdata( 'nci_60.txt' );

data = rawdata.data;

cancer_factor = cell( length( rawdata.textdata ), 1 );
for i=1:length( rawdata.textdata )

    cancer_factor{i} = rawdata.textdata{i}(2:4);
    
end

% running SVD/PCA
X = data - repmat( mean( data ), size( data, 1 ), 1 );

[ U S V ] = svd( X, 0 );

% ploting samples based on two first components, without any coloring.
fig1 = figure(1);
plot( V(:,1), V(:,2), '.', 'MarkerSize', 6, 'LineWidth', 3  )
hold on
text( V(:,1), V(:,2) + 0.01, cancer_factor, 'FontSize', 8 )
hold off
xlabel( 'pc1' )
ylabel( 'pc2' )
% saveas(fig1, '../report/images/q34pcs1', 'epsc')

% hierarchical clustering
% D = pdist( data', 'euclidean' );
% z = linkage( D, 'average' );
% fighc = figure(2020);
% dendrogram( z, 'labels', cancer_factor )
% dendrogram( z, 'orientation', 'left' )
% idx = str2num( get( gca, 'YTickLabel' ) ); %#ok
% set( gca, 'YTickLabel', cancer_factor( idx ) )
% saveas(fighc, '../report/images/q35hc', 'epsc')

% two-way clustering
% var_thres = 4;
% data_filtered = data(var( data, 0, 2 ) > var_thres,:);
% fighm = figure(2101);
% clustergram( data_filtered, 'ColumnLabels', cancer_factor )

% kmeans clustering, seeds to make them the same.
rand( 'seed', 0 ); %#ok
randn( 'seed', 0 ); %#ok

% Makes a 2d plot based on the first two principal components. Running
% kmeans clustering on it using k = 9
% k = 9;
% cl = kmeans( data', k );

% + 0.015 is just to position the text a bit above the points just plotted,
% because of the additional LineWidth parameter added to scatter.
% fig2 = figure(2);
% scatter( V(:,1), V(:,2), 30, cl, 'LineWidth', 3 )
% hold on
% text( V(:,1), V(:,2) + 0.015, cancer_factor, 'FontSize', 8 )
% hold off
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% title([ 'Coloring using k = ', num2str(k) ]);
% saveas(fig2, '../report/images/q34pcs2', 'epsc')

% Makes a 2d plot based on the first two principal components. Running
% kmeans clustering on it using k = 6.
% k = 6;
% cl = kmeans( data', k );

% + 0.015 is just to position the text a bit above the points just plotted,
% because of the additional LineWidth parameter added to scatter.
% fig3 = figure(3);
% scatter( V(:,1), V(:,2), 30, cl, 'LineWidth', 3  )
% hold on
% text( V(:,1), V(:,2) + 0.015, cancer_factor, 'FontSize', 8 )
% hold off
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% title( [ 'Coloring using k = ', num2str(k) ]);
% saveas(fig3, '../report/images/q34pcs3', 'epsc')

% This makes a 3d plot based on the three first principal components.
% Running kmeans clustering on it using k = 7
% k = 7;
% cl = kmeans( data', k );

% fig47 = figure(47);
% scatter3( V(:,1), V(:,2), V(:,3), 30, cl, 'LineWidth', 3  )
% hold on
% text( V(:,1), V(:,2) - 0.01, V(:,3), cancer_factor, 'FontSize', 8 )
% hold off
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% zlabel( 'pc3' )
% title([ 'Coloring using k = ', num2str(k) ]);
% saveas(fig47, '../report/images/q34pcs2', 'epsc')

% Makes a 2d plot based on the first two principal components. Running
% kmeans clustering on it using the vector k = [1,...,9]
% k = 2:10;
% 
% for i = 1:length(k)
%     subfigkmeans = subplot(3,3, i);
%     
%     cl = kmeans( data', k(i) );
%     
%     % + 0.015 is just to position the text a bit above the points just plotted,
%     % because of the additional LineWidth parameter added to scatter.
%     scatter( V(:,1), V(:,2), 30, cl, 'LineWidth', 2 )
%     hold on
%     text( V(:,1), V(:,2) + 0.015, cancer_factor, 'FontSize', 8 )
%     hold off
%     xlabel( 'pc1' )
%     ylabel( 'pc2' )
%     title([ 'Coloring using k = ', num2str(k(i)) ]);
% end
% saveas(subfigkmeans, '../report/images/q36pcs5', 'epsc')

% cl <- kmeans(t(data), 6)
% ploting of the samples, plotting points type is determined by clusters, color is determined by cancer type
% plot(V[,1], V[,2], pch = cl$cluster) 
% % adding the labels
% text ( V[,1], V[,2]+0.01, colnames(data), cex = 0.55, col = cancer_colors)
 
% density estimation
 
%  histogram select the number of bins, as long as it sums to nine in
%  total.
bins_vector = 3:3:27;

% Question 3.7a
% figure(4),
% for i=1:length(bins_vector)
%     [ xout n ] = Normlized_Histogram(V(:, 2), bins_vector(i));
%     subfig1 = subplot(3, 3, i);
%     bar(xout, n, 'BarWidth', 1);
%     title(['Number of bins = ' int2str(bins_vector(i))])
% end

% saveas(subfig1, '../report/images/q37histograms', 'epsc')

% figure(4)
% hist( V(:,2), bins )
% xlabel( 'pc2' )

%  kernel density estimator
%  let bw unespecified to select it automatically
bw = 0.02:0.01:0.1;  % bw = Bin Width

% Question 3.7b
% figure(5),
% for i=1:length(bw)
%     subfig2 = subplot(3,3, i);
%     ksdensity( V(:,2), 'width', bw(i))
%     title(['Bandwidth = ' num2str(bw(i))]) 
% end
% saveas(subfig2, '../report/images/q37kde', 'epsc')

% Let ksdensity find it own optimal bandwidth
% [~,~,u] = ksdensity( V(:,2));
% fig15 = figure(15),
% ksdensity( V(:,2) )
% title(['Bandwidth = ' num2str(u)])


% kernel density 2d
% the bandwidth is selected automatically
% n is the number of points used to make the two-dimensional grid
n = 100;
% The bandwidth we run through, we overshoot by 2, which does not get
% plottet.
bandwidth = [0.001:0.001:0.01;0.001:0.001:0.01];

% [ h f X Y ] = kde2d( V(:,1:2), n );
% figure(155)
% for i=1:4
%     [ h f X Y ] = kde2d( V(:,1:2), n, bandwidth(1,i), bandwidth(2,i));
%     subfig3 = subplot(2,2, i);
%     surf( X, Y, f )
%     xlabel( 'pc1' )
%     ylabel( 'pc2' )
%     zlabel( 'density' )
%     title(['Bandwidth = ' num2str(bandwidth(1,i))]) 
%     axis tight
% end
% saveas(subfig3, '../report/images/q373dkde1', 'epsc')
% 
% figure(156)
% for i=1:4
%     [ h f X Y ] = kde2d( V(:,1:2), n, bandwidth(1,i+4), bandwidth(2,i+4));
%     subfig4 = subplot(2,2, i);
%     surf( X, Y, f )
%     xlabel( 'pc1' )
%     ylabel( 'pc2' )
%     zlabel( 'density' )
%     title(['Bandwidth = ' num2str(bandwidth(1,i+4))]) 
%     axis tight
% end
% saveas(subfig4, '../report/images/q373dkde2', 'epsc')

% [ h f X Y ] = kde2d( V(:,1:2), n );

% 3d plot if wanted
% figure(1000),
% surf( X, Y, f )
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% zlabel( 'density' )
% axis tight

% figure(6)
% % par is the number of contours to plot
% par = 10;
% contour( X, Y, f, par )
% hold on
% plot( V(:,1), V(:,2), 'rx' )
% hold off
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% axis tight
% grid on

%  mixtures of gaussians
%  G is the number of components, default G = 3
G = 2:1:5;

% Option set for guassian mixture model class.
options = statset( 'Display', 'off', 'MaxIter', 500 );

% grid size for the meshgrid function.
gsize = 100;

% Makes a 3d subplot based on G from above, these are the ones that gives a
% pretty smooth density without extreme spikes.
% for i=1:length(G)
%     omog = gmdistribution.fit( V(:,1:2), G(i), 'Options', options, 'Replicates', 3 );
%     
%     [ X Y ] = meshgrid( linspace( -0.4, 0.4, gsize ), linspace( -0.4, 0.4, gsize ) );
%     mf = pdf( omog, [ X(:) Y(:) ] );
%     
%     figure(200)
%     subfig5 = subplot(2,2,i);
%     surf( X, Y, reshape( mf, gsize, gsize ) )
%     xlabel( 'pc1' )
%     ylabel( 'pc2' )
%     zlabel( 'density' )
%     title(['Number of components = ' num2str(G(i))])
%     axis tight
% end
% saveas(subfig5, '../report/images/q37mog', 'epsc')

% figure(7)
% contour( X, Y, reshape( mf, gsize, gsize ), par )
% hold on
% plot( V(:,1), V(:,2), 'rx' )
% hold off
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% axis tight
% grid on

% to select G automatically loop over options for G ans select the one with
% better (minimum) score (aic)
G = 1:12;

aic = zeros( length( G ), 1 );

for i=1:length( G )
    omog = gmdistribution.fit( V(:,1:2), G(i), 'Options', options, 'Replicates', 3 );
    
    aic(i) = omog.AIC;
end

[ baic bestG ] = min( aic );
fprintf('The best number of gaussian components in the mixture model is: %d, this is with baic %f\n', bestG, baic)

omog = gmdistribution.fit( V(:,1:2), bestG, 'Options', options, 'Replicates', 3 );

[ X Y ] = meshgrid( linspace( -0.4, 0.4, gsize ), linspace( -0.4, 0.4, gsize ) );
mf = pdf( omog, [ X(:) Y(:) ] );

% fig210 = figure(210);
% surf( X, Y, reshape( mf, gsize, gsize ) )
% xlabel( 'pc1' )
% ylabel( 'pc2' )
% zlabel( 'density' )
% title(['Number of components = ' num2str(bestG)])
% axis tight
% saveas(fig210, '../report/images/q37mogbest', 'epsc')