function plot_results( prob_map, centroid, covariance, points )
% plot_contour(prob_im,centroid,covariance): plot the centroid and the
% contours of the corresponding gaussian on top of the model
% probabability map
%
%
% centroid, and the contours of the probability map (prob_map)
% prob_map: a M x N data matrix
% centroid: 2 x 1 vector
% covariance: 2 x 2 matrix

if nargin < 4
	points = 200;
end

dx = linspace( 1, size( prob_map, 1 ), points );
dy = linspace( 1, size( prob_map, 2 ), points );

[ grdx grdy ] = meshgrid( dx, dy );
dens = mvnpdf( [ grdx(:) grdy(:) ], centroid', covariance );
dens = reshape( dens, length( dx ), length( dy ) );

figure
imagesc( prob_map )
hold on
plot( centroid(1), centroid(2), 'rx', 'MarkerSize', 100, 'Linewidth', 2 )
contour( dx, dy, dens, 'g-' )
hold off
axis tight
xlabel( 'x' )
ylabel( 'y' )

%im = zeros(size(prob_map));
%Z = 1/( (2*pi) * det(covariance)^(1/2));

%for i = 1:size(im,1)
%  for j = 1:size(im,2)
%    x = [i; j];
%    im(i,j) = Z * exp(-1/2*(x-centroid)'*inv(covariance)*(x-centroid));
%  end
%end

%figure()
%imagesc(prob_map);
%hold on;
%plot(centroid(2),centroid(1),'rx','MarkerSize',100,'LineWidth',2)
%contour(im,'g-');
