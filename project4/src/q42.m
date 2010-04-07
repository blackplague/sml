function q42( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    function [m1 m2 C1 C2 train_err] = q42trainclass( path )
        
        % As we consider a two-class problem we need to model the
        % probabilistic density 
        % p(C1|x) = p(x|C1)p(C1) / sum_i(p(x|Ci)p(Ci))

        addpath( './hfunc/' ) % Helper functions
        
        [ userfeatures rv ] = getuserinfo( path );
        [uf dr] = discretize( userfeatures, rv );

        m1 = mean(uf((dr == 1), 3:end));
        m2 = mean(uf((dr == 2), 3:end));
        
        C1 = cov(uf((dr == 1), 3:end));
        C2 = cov(uf((dr == 2), 3:end));
        
        p = zeros(length(dr), 1);
        
        % Calculate non-normalized probability of the two classes given
        % data
        for i = 1:length(dr)
                p(i,1) = mvnpdf( uf(i,3:end), m1, C1 );
                p(i,2) = mvnpdf( uf(i,3:end), m2, C2 );
        end
       
        % The hacked way of getting the normalization coefficient for the
        % distributions.
        normcoeff1 = sum(p(:,1)); normcoeff2 = sum(p(:,2));
        p(:,1) = p(:,1)/normcoeff1; p(:,2) = p(:,2)/normcoeff2;
        
        pdr = predicteddiscreterating(p);
        
        train_err = misclassificationrate(dr, pdr);
        
    end

    function err_rate = misclassificationrate(dr, pdr)

        misclass = 0;        
        
        % Loop to find misclassified
        for i=1:length(dr)
            if(dr(i) ~= pdr(i))
                misclass = misclass + 1;
            end
        end

        err_rate = (misclass / length(dr))*100;        
    end

    function pdr = predicteddiscreterating(p)
        % Init all to something that not a class
        pdr = 2*ones(length(p),1);
        
        % Assign class by using the max of probability of the two classes
        for i=1:length(p)
            [~, I] = max(p(i,:));
            pdr(i) = I;
        end
    end

    function test_err = classifier(m1, m2, C1, C2, path, Theta)
        
        [ userfeatures rv ] = getuserinfo( path );        
        [uf dr] = discretize( userfeatures, rv );
                
        p = zeros(length(dr), 2);
        
        for i = 1:length(dr)
                p(i,1) = mvnpdf( uf(i,3:end), m1, C1 );
                p(i,2) = mvnpdf( uf(i,3:end), m2, C2 );
        end

        % The hacked way of getting the normalization coefficient for the
        % distributions.
        normcoeff1 = sum(p(:,1)); normcoeff2 = sum(p(:,2));
        p(:,1) = p(:,1)/normcoeff1; p(:,2) = p(:,2)/normcoeff2;
        
        pdr = predicteddiscreterating(p);
        
        test_err = misclassificationrate(dr, pdr);
        
    end

    % Function for discretizing the rating matrix, it maps 1 and 2 -> 0
    % for 'dislike', 4 and 5 -> 1 for 'like' and just throws the rest out.
    function [newUF newDR] = discretize( uf, rv )
        
        N = length(find(rv == 1 | rv == 2 | rv == 4 | rv == 5));
        
        newUF = zeros(N, size(uf,2));
        newDR = zeros(N,1);
        
        c = 1;
        for i=1:length(rv)
            switch rv(i,1)
                case {1,2}
                    newDR(c,1) = 2;
                    newUF(c,:) = uf(i,:);
                    c = c + 1;
                case {4,5}
                    newDR(c,1) = 1;
                    newUF(c,:) = uf(i,:);
                    c = c + 1;
            end
        end
        
    end

    function run()
        
        % Return vector for values from the train classifier function
        % ret_vec_tr <==> [m1 m2 C1 C2 tr_err]
        m1 = zeros(5,24); m2 = zeros(5,24); C1 = cell(5,1); C2 = cell(5,1);
        tr_err = zeros(5,1);
        
        for i=1:size(tr_err,1)
            [m1(i,:) m2(i,:) C1{i} C2{i} tr_err(i)] = q42trainclass(['u' num2str(i) '.base']);
            fprintf('Done training %d ', i)
        end

        fprintf('\n')
        for i=1:size(tr_err,1)
            fprintf('Train error %d: %2.5f\n', i, tr_err(i,1))
        end
        
        te_err = zeros(5,1);
        
        for i=1:size(te_err,1)
            te_err(i) = classifier(m1(i,:), m2(i,:), C1{i}, C2{i}, ...
                ['u' num2str(i) '.test']);
            fprintf('Done classifying %d ', i)
        end
        
        fprintf('\n')
        
        for i=1:size(te_err,1)
            fprintf('Test error %d: %2.5f\n', i, te_err(i,1))
        end        
        
    end
    run();
end