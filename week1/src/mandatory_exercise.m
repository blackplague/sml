function [ output_args ] = mandatory_exercise( input_args )
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    function [ y ] = exercise11(mu, sigma, number)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
        
        x = randn(2, number);

        L = chol(sigma, 'lower');
            
        y = (mu + L * x)';
        
    end

    function [ mu_ml ] = exercise12a( x )
        %UNTITLED3 Summary of this function goes here
        %   Detailed explanation goes here
        
        mu_ml = sum(x',2)/length(x);

    end

    function [ sigma_est ] = exercise12b( mu_ml, y, number )
        %UNTITLED4 Summary of this function goes here
        %   Detailed explanation goes here

        y = y';

        N = length(y);
    
        numberOnes = ones(number,2);
    
        tmp = [mu_ml(1) * numberOnes(:,1)'; mu_ml(2) * numberOnes(:,2)'];
        
        sigma_est = 1 / (N - 1) * (y - tmp)*(y - tmp)';

    end

    function [ output_args ] = exercise12c( correct_mean, estimated_mean, data_points )
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here

        figure,
        hold on
        scatter(data_points(:,1), data_points(:,2), 'xg')
        scatter(correct_mean(1), correct_mean(2), 'o', 'MarkerFaceColor', 'b')
        scatter(estimated_mean(1), estimated_mean(2), 'o', 'MarkerFaceColor', 'r')
        hold off
        
%         figure,
%         hold on
%         plot(data_points(:,1),data_points(:,2), 'o', 'MarkerFaceColor', 'g')
%         plot(correct_mean(1), correct_mean(2), 'o', 'MarkerFaceColor', 'r')
%         plot(estimated_mean(1), estimated_mean(2), 'o', 'MarkerFaceColor', 'b')
%         hold off
        
    end

    function [ output_args ] = exercise13( y, bins )
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here

        x1 = y(:,1);
        x2 = y(:,2);

        figure,
        hold on
        subplot(1,2,1)
        hist(x1, bins)
        subplot(1,2,2)
        hist(x2, bins)
        hold off
        
    end

    function [ output_args ] = exercise14( y, sigma )
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here

        % Numbers of new datapoints
        N1 = 100;
        N2 = 1000;
        N3 = 10000;
        
        % Generate the new datapoints
        y1 = exercise11(N1, sigma, N1);
        y2 = exercise11(N2, sigma, N1);
        y3 = exercise11(N3, sigma, N1);        
        
        % Plot the original data.
        figure, hist3( y )
        
        figure, hist3( y1 )
        figure, hist3( y2 )
        figure, hist3( y3, [10 10] )
        figure, hist3( y3, [15 15] )
        figure, hist3( y3, [20 20] )
    
    end

    function exercise15( samples )
       
        results = zeros(3, 1);
        
        for i=1:length( samples )
            results(i) = calculator15( samples(i) );
        end
        
        figure, 
        hold on
        plot(samples, results, '--x')
        hold off
        
    end

    function [ diff_mu ] = calculator15( samples )
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here

        fprintf(1, 'Number of samples: %i\n', samples)
        
        repetitions = 1000;

%         mu_x = zeros(repetitions, samples);
%         y_hat = zeros(repetitions, samples);
        
        mu_x = 0; y_hat=0;

        for i=1:repetitions
            x = generateSamples( samples );
            mu_x = mu_x + mean(x);
            y = transform(x(:), 1/2);
            y_hat = y_hat + (sum(y)/samples);
        end
        
        mean_x = mu_x/repetitions;
        mean_y = y_hat/repetitions;

        fprintf(1,'Mean of x: %f\n', mean_x)
        fprintf(1,'Mean of y: %f\n', mean_y)
        
        mu_y = 1/mean_x;
        
        diff_mu = abs(mu_y - mean_y);
        
        fprintf(1, 'Difference between mu_y and mean_y: %f\n\n', diff_mu)

    end

        function [ x ] = generateSamples( samples )
            
            x = rand(1, samples);
            
        end

        function [ y ] = transform( z, lambda )
            y = -lambda^(-1)*log(1-z);
        end

    function exercise16( Sigma, samplessize, samplertype, Stepsize )
        
        % case 1 Statistical Methods for Machine Learning
        % MCMC for multivariate normal
        
        mu = [ 1 ; 1 ];
        
        
        sampler = samplertype; % options: {'exact','MH'}
        L = samplessize; % number of samples
        stepsize = Stepsize; % step-size in M-H proposal density
        lagmax = 50; % maximum lag in autocorrelation
        
        y = zeros(length(mu),L); % samples
        C = chol(Sigma,'lower');
        
        switch sampler
            case 'exact'
                
                y = mu(:,ones(L,1)) + C * randn(length(mu),L);
                
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
                            + 0.5*(y(:,l-1)-mu)'*Sigmainv*(y(:,l-1)-mu))
                        
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
        
        R = zeros(length(mu),lagmax);
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
    end
    
    function exercise18( image_path )
       
        original_image = im2double (imread( image_path ));
        
        t_data = extract_training_region( original_image );
        
        ts = size(t_data);
        
        % Calculate mean vector for RGB in test region
        mu_ml = (1 / length(t_data)) * sum(t_data);
        numberOnes = ones(ts(1), ts(2));
    
        tmp = [mu_ml(1) * numberOnes(:,1)'; mu_ml(2) * numberOnes(:,2)' ...
            ; mu_ml(3) * numberOnes(:,3)'];
        
        sigma = cov(t_data);

        %RGB = im2double (imread ( './images/kande1.jpg' ));

        new_image = getProbability(original_image, mu_ml, sigma);
        
        figure, imshow(new_image);
                
    end

    function [ new_image ] = getProbability( image, mu_ml, sigma)

        sqr_det_sigma = sqrt(det(sigma));
        
        inverse_sigma=inv(sigma);
        
        [sm sn u] = size(image);
        
        new_image = zeros(sm,sn);
        
        % Needs normalization
        for m=1:sm
            for n=1:sn
                x = [image(m,n,1); image(m,n,2); image(m,n,3)];
                
                new_image(m,n) = 1/(2*pi)^(3/2)*1/sqr_det_sigma * ...
                    exp(-0.5*(x - mu_ml')'*inverse_sigma*(x - mu_ml'));
            end
        end
    end

    function [ training_data ] = extract_training_region( image )
       
        %% Read Image
        
        %% Define training and test regions
        tll = [328, 150]; % training: lower left
        tur = [264, 330]; % training: upper right
        
        %% Show image and regions
        figure,
        hold on
        plot ([tur(2), tll(2), tll(2), tur(2), tur(2)], ...
            [tur(1), tur(1), tll(1), tll(1), tur(1)], ...
            'k-', 'linewidth', 2);
        text (tll(2), tur(1)-15, 'Training');
        
        hold off
        
        %% Extract regions
        training_data = image(tur(1):tll(1), tll(2):tur(2), :);
        training_data = reshape (training_data, size (training_data, 1) * size (training_data, 2), 3);
        
    end


    function run( ~ )
        
%         numberOfPoints = 100;
        
%         mu = 1.0;
        
%         sigma = [ 0.3 0.2; 0.2 0.2 ];
        
%         y = exercise11( mu, sigma, numberOfPoints );
        
%         [ mu_ml ] = exercise12a( y );
        
%         sigma = exercise12b( mu_ml, y, numberOfPoints );
    
%         exercise12c( [mu mu], mu_ml, y );
        
%         exercise13( y, 5 );
        
%         exercise14( y, sigma )
        
%         samples = [ power(10,1); power(10,2); power(10,3) ];
        
%         exercise15( samples )
                
        % samplertype is either exact or MH, stepsize = 0.5
%         samplessize = 500; samplertype = 'exact'; stepsize = 0.7;
%         exercise16( sigma, samplessize, samplertype, stepsize )
        
        image_path = './images/';
        image_name = 'kande1.jpg';

        exercise18( strcat( image_path, image_name ))
        
    end

    run()

end