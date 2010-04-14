function [ output_args ] = q43( input_args )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    function [ RMSE Y ] = mf( X, K, lambda )
         
        for i = 1:length(K)
            for j = 1:length(lambda)
                [ Y, U, V, RMSE ] = matrixfactorization(X, K(i), lambda(j));
                fprintf('Training, K=%d, lambda=%d, RMSE=%2.4f\n', K(i), lambda(j), RMSE)
            end            
        end
    end

    function [ RMSEte ] = test( X, Y )
        
        indx = X > 0 ;
        Ntr = sum(sum(indx));
        
        E = X - Y ; % compute residual for all entries
        
        RMSEte = sqrt(sum(sum((indx .* E).^2))/Ntr);
        
    end

    function run()

        K = [5 10 25 50 100]; lambda = [0 0.1 0.4 1 2 4 8 16 25];
        
        % Colors for plotting the RMSE for training, test and cross
        % validation.
        c = [ ...
            0 0 0; 1 0 0; 0 1 0; 0 1 1; 0 0 1; ...
            1 1 0; 1 0 1; 1 0.3 0.6; 0.3 1 0.6; 0.7 1 0.2 ...
            ];
        
        file_exists = zeros(5,2);
        
        for i = 1:length(file_exists)
           file_exists(i,1) = exist(['./data/rmse_tr' num2str(i) '.mat' ], 'file' );
           file_exists(i,2) = exist(['./data/rmse_te' num2str(i) '.mat' ], 'file' );
        end
        
        RMSEmean_tr = zeros(length(K),length(lambda));
        RMSEmean_te = zeros(length(K),length(lambda));
        
        if(isempty(find(~file_exists,1)))
            
            for i = 1:length(file_exists)
                load(['./data/rmse_tr' num2str(i) '.mat' ])
                load(['./data/rmse_te' num2str(i) '.mat' ])
                RMSEmean_tr = RMSEmean_tr + RMSE;
                RMSEmean_te = RMSEmean_te + RMSEte;
                
                fig = figure(i);
                hold on
                for j=1:length(K)
                    plot(lambda(:), RMSE(j,:), '-', 'LineWidth', 1, 'Color', c(j,:))
                    plot(lambda(:), RMSEte(j,:), '-', 'LineWidth', 1, 'Color', c(length(K)+j, :))                    
                end
                legend( ...
                    ['K = ' num2str(K(1)) ' train'], ...
                    ['K = ' num2str(K(1)) ' test'],  ...
                    ['K = ' num2str(K(2)) ' train'], ...
                    ['K = ' num2str(K(2)) ' test'],  ...
                    ['K = ' num2str(K(3)) ' train'], ...
                    ['K = ' num2str(K(3)) ' test'],  ...
                    ['K = ' num2str(K(4)) ' train'], ...
                    ['K = ' num2str(K(4)) ' test'],  ...
                    ['K = ' num2str(K(5)) ' train'], ...
                    ['K = ' num2str(K(5)) ' test'], ...
                    'Location', 'SouthEast')
                title(['RMSE training and test set ' num2str(i)])
                xlabel('Lambda')
                ylabel('RMSE')
                hold off
            saveas(fig, ['../report/images/q43RMSEtrte' num2str(i)], 'epsc')
            end
            
            RMSEmean_tr = RMSEmean_tr/length(K);
            RMSEmean_te = RMSEmean_te/length(K);
            
            fig6 = figure(6);
            hold on
            for i=1:length(K)
                plot(lambda(:), RMSEmean_tr(i,:), '-', 'LineWidth', 1, 'Color', c(i,:))
                plot(lambda(:), RMSEmean_te(i,:), '-', 'LineWidth', 1, 'Color', c(length(K)+i, :))
            end
            legend( ...
                ['K = ' num2str(K(1)) ' train'], ...
                ['K = ' num2str(K(1)) ' test'],  ...
                ['K = ' num2str(K(2)) ' train'], ...
                ['K = ' num2str(K(2)) ' test'],  ...
                ['K = ' num2str(K(3)) ' train'], ...
                ['K = ' num2str(K(3)) ' test'],  ...
                ['K = ' num2str(K(4)) ' train'], ...
                ['K = ' num2str(K(4)) ' test'],  ...
                ['K = ' num2str(K(5)) ' train'], ...
                ['K = ' num2str(K(5)) ' test'], ...
                'Location', 'SouthEast')
            title(['RMSE Cross validation for training and test'])
            xlabel('Lambda')
            ylabel('RMSE')
            hold off
            saveas(fig6, '../report/images/q43RMSEmean', 'epsc')
        end
        
%         if(exist('mfvar.mat', 'file') ~= 2)
%             fprintf('mfvar.mat not found performing matrixfactorization to save variables Y_mf U and V.')
%             [ Xt, timestamps, user, movie, genre ] = ...
%                 readmovielens('ml-data/', [ 'u' '.data']);
%             K=50;lambda=5;
%             [Y_mf U V ~] = matrixfactorization(Xt,K,lambda);
%             save('mfvar.mat', 'Y_mf', 'U', 'V');
%         else
%             load('mfvar.mat');
%             fprintf('mfvar.mat found loaded variables Y_mf, U and V.')
%         end
%         
%         K = [5 10 25 50 100]; lambda = [0 0.1 0.4 1 2 4 8 16 25];
%         Xtr = cell(5,1);
%         Xte = cell(5,1);
%         
%         % Read all ratings matrices into Xr cell.
%         for i=1:length(Xtr)
%             [ Xt, timestamps, user, movie, genre ] = ...
%                 readmovielens('ml-data/', [ 'u' num2str(i) '.base']);
%             Xtr{i} = Xt;%(1:100,1:100);
%             [ Xtest, timestamps, user, movie, genre ] = ...
%                 readmovielens('ml-data/', [ 'u' num2str(i) '.test']);
%             Xte{i} = Xtest;%(1:100,1:100);            
%         end
%         
% %         RMSEcross_tr = zeros(length(K), length(lambda));
% %         RMSEcross_te = zeros(length(K), length(lambda));
% 
%         RMSE = zeros(length(K), length(lambda));
%         RMSEte = zeros(length(K), length(lambda));
% 
% %         save(['/vol/tmp/bp/mfvar_tr' num2str(i) '_k' num2str(K(m)) '_l' ...
% %             num2str(lambda(n)) '.mat'], 'Y', 'U', 'V')
%         
%         for i=1:length(Xtr)
%             for m=1:length(K)
%                 for n=1:length(lambda)
%                     [ Y U V RMSE(m,n) ] = matrixfactorization(Xtr{i}, ...
%                         K(m), lambda(n));
%                     save(['/vol/tmp/bp/mfvar_tr' num2str(i) '_k' ...
%                         num2str(K(m)) '_l' num2str(lambda(n)) ...
%                         '.mat'], 'Y', 'U', 'V')
%                     fprintf('Done factorizing and saving for training %d, K = %d, lambda = %d, RMSE = %d\n', i, K(m), lambda(n), RMSE(m,n))
%                     RMSEte(m,n) = test(Xte{i}, Y);
%                 end
%             end
%             save(['/vol/tmp/bp/rmse_tr' num2str(i) '.mat'], 'RMSE')
%             save(['/vol/tmp/bp/rmse_te' num2str(i) '.mat'], 'RMSEte')
%             fprintf('Done with RMSE for training and test (%d)\n', i)
%         end
%         
% 
%         
% %         for i=1:length(Xtr)
% %             
% %             Y = cell(length(K), length(lambda));
% %             
% % %             figure(i)
% % %             hold on
% %             
% %             for m=1:length(K)
% %                 for n=1:length(lambda)
% %                     [ RMSE(m,n) Y{m,n} ] = mf(Xtr{i}, K(m), lambda(n));
% %                     RMSEte(m,n) = test(Xte{i}, Y{m,n});
% %                 end
% %             end
%             
% %             RMSEcross_tr = RMSEcross_tr+RMSE;
% %             RMSEcross_te = RMSEcross_te+RMSEte;
% %             
% %             plot(lambda(:), RMSE(1,:), '-', 'LineWidth', 1, 'Color', 'blue')
% %             plot(lambda(:), RMSE(2,:), '-', 'LineWidth', 1, 'Color', 'red')
% %             plot(lambda(:), RMSEte(1,:), '-', 'LineWidth', 1, 'Color', 'yellow')
% %             plot(lambda(:), RMSEte(2,:), '-', 'LineWidth', 1, 'Color', 'black')
% %             legend(['K = ' num2str(K(1)) ' train'], ['K = ' ... 
% %                 num2str(K(2)) ' train'], ['K = ' num2str(K(1)) ...
% %                 ' test'], ['K = ' num2str(K(2)) ' test'])
% %             title(['RMSE Training set ' num2str(i)])
% %             xlabel('Lambda')
% %             ylabel('RMSE')
% %             hold off
% %             saveas(i, ['../report/images/q43_' num2str(i)], 'epsc')
% %         end
%         
% %         figure(100)
% %         hold on
% %         plot(lambda(:), RMSEcross_tr(1,:)/5, '-', 'LineWidth', 1, 'Color', 'blue')
% %         plot(lambda(:), RMSEcross_tr(2,:)/5, '-', 'LineWidth', 1, 'Color', 'red')
% %         plot(lambda(:), RMSEcross_te(1,:)/5, '-', 'LineWidth', 1, 'Color', 'yellow')
% %         plot(lambda(:), RMSEcross_te(2,:)/5, '-', 'LineWidth', 1, 'Color', 'black')
% %         legend(['K = ' num2str(K(1)) ' train'], ['K = ' ...
% %             num2str(K(2)) ' train'], ['K = ' num2str(K(1)) ...
% %             ' test'], ['K = ' num2str(K(2)) ' test'])
% %         title('Cross validation - avg. RMSE for training and test')
% %         xlabel('Lambda')
% %         ylabel('RMSE')
% %         hold off
% %             saveas(100, ['../report/images/q43_' num2str(100)], 'epsc')
    end
    run();
end