function q41( ~ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function [RMSE w_ml] = q41train( path )

        addpath('./hfunc/')
        
        [uf, rv, ~] = getuserinfo( path );
        
%         normage = normalizedata(uf(:,3));
%         normtime = normalizedata(uf(:,8));
        
        phifn=@(x)[1 x];

%         Phi1 = generate_designmatrix([ normage uf(:,4:end-1) normtime], phifn);
%         Phi1 = generate_designmatrix([ normage uf(:,4:end)], phifn);
        Phi1 = generate_designmatrix(uf(:,3:end), phifn);
        
        w_ml = inv(Phi1'*Phi1)*Phi1'*rv;
        
        y_1 = Phi1 * w_ml;

        % Truncate values over 5 to 5, and below 0 to 0.
        y_1(y_1 > 5) = 5; y_1(y_1 < 0) = 0;
        
        RMSE = sqrt(sum((rv-y_1).^2)/length(y_1));
        
    end

    function [RMSE] = q41test(w_ml, path)
        
        [uf, rv, ~] = getuserinfo( path );

        phifn=@(x)[1 x];
        Phi1 = generate_designmatrix(uf(:,3:end), phifn);        
        
        y_1test = Phi1 * w_ml;
        
        y_1test(y_1test > 5) = 5; y_1test(y_1test < 0) = 0;
        
        RMSE = sqrt(sum((rv-y_1test).^2)/length(y_1test));
        
    end

    function [RMSE] = q41testmean(mean_vector, rm)
        
        mean_matrix = repmat(mean_vector,size(rm,1),1);
        
        indx = rm>0;
        Ntr = sum(sum(indx));
        
        RMSE = sqrt(sum(sum(((indx .* rm) - (indx .* mean_matrix)).^2))/Ntr);
    end

    function normalizedvector = normalizedata( vector )
    %normalizedata recieves a vector of data and return a normalized vector
    %by performing the following normalization data((i) - mean(data)) /
    %std(data) (Whitening of the data).
        
        meanofdata = mean(vector);
        stdofdata = std(vector);
        
        normalizedvector = ...
            (vector - repmat(meanofdata,size(vector,2),1)) * (1/stdofdata);
    
    end
        
    function run()
        
        RMSEtrain = zeros(5,1);
        RMSEtest = zeros(5,1);
        w_ml = zeros(25,5);

        for i=1:length(RMSEtrain)
            fprintf('Training %d\n', i)
            [ RMSEtrain(i,1) w_ml(:,i) ] = q41train( [ 'u' num2str(i) '.base'] );
            [ RMSEtest(i,1) ] = q41test(w_ml(:,i), ['u' num2str(i) '.test']);
            fprintf('Done\n')
        end
        
        fprintf('RMSE for training:\n')
        for i=1:length(RMSEtrain)
            fprintf('RMSE%d: %2.6f\t ', i, RMSEtrain(i,1) )
        end
        fprintf('\n')

        fprintf('RMSE for testing:\n')
        for i=1:length(RMSEtrain)
            fprintf('RMSE%d: %2.6f\t ', i, RMSEtest(i,1) )
        end
        fprintf('\n')
        
        fprintf('Mean RMSE training: %2.6f\n', mean(RMSEtrain))
        fprintf('Mean RMSE test: %2.6f\n', mean(RMSEtest))        
        
        for i=1:length(w_ml)
            fprintf('w%d\t', i)
        end
        fprintf('\n')
        for i=1:size(w_ml,2)
            for j=1:length(w_ml)
                fprintf('%2.4f\t', w_ml(j,i))
            end
            fprintf('\n')
        end
        
        RMSEtest_avg = zeros(5,1);
        w_ml_mean = mean(w_ml,2);
        for i = 1:length(RMSEtest_avg)
            [ RMSEtest_avg(i,1) ] = q41test(w_ml_mean, ['u' num2str(i) '.test']);
        end

        fprintf('RMSE for testing based on avg. w_ml:\n')
        for i=1:length(RMSEtest_avg)
            fprintf('RMSE%d (avg): %2.6f\t ', i, RMSEtest_avg(i,1) )
        end
        fprintf('\n')

        fprintf('Mean RMSE test (Based on mean w_ml): %2.6f\n', mean(RMSEtest_avg))        
        
        mean_vec = zeros(5,1682);
        RMSE_mean = zeros(5,1);
        for i=1:5
            [~, ~, rm_train] = getuserinfo( ['u' num2str(i) '.base'] );
            [~, ~, rm_test ] = getuserinfo( ['u' num2str(i) '.test'] );
            
            mean_vec(i,:) = meannonzero( rm_train );
            RMSE_mean(i) = q41testmean(mean_vec(i,:), rm_test);
        end
        
        for i=1:length(RMSE_mean)
            fprintf('RMSE%d based on average movie rating: %2.10f\n', i, RMSE_mean(i))
        end
        
        fprintf('Mean RMSE: %2.10f\n', mean(RMSE_mean))
    end
run();
end