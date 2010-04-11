function [ output_args ] = q46( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    function q46knn(k,space,textOn)
    %q46knn Performs knn on either U or V, the space used depends on the
    %variable space, which can be either 'U' or 'V'. By default the space
    %used is U and if no k is given the default value is 20.
                
        if nargin < 1
            k = 20;
        end            
        
        if nargin < 2
            space = 'U';
        end
        
        if nargin < 3
            textOn = 'false';
        end
        
        if(exist('mfvar.mat', 'file') ~= 2)
            fprintf('mfvar.mat not found performing matrixfactorization to save variables Y_mf U and V.')
            [ Xt, timestamps, user, movie, genre ] = ...
                readmovielens('ml-data/', [ 'u' '.data']);
            K=50;lambda=5;
            [Y_mf U V ~] = matrixfactorization(Xt,K,lambda);
            save('mfvar.mat', 'Y_mf', 'U', 'V');
        else
            load('mfvar.mat');
            fprintf('mfvar.mat found loaded variables Y_mf, U and V.')
        end
        
        U_old = U - repmat( mean( U ), size( U, 1 ), 1 );
        
        [U_new S V] = svd(U_old, 0);
        
        path = 'u.data';
        
        % Load once, use in all index functions, alternative, send specific
        % path to index function to get indices in training or test sets.
        [~, ~, user, ~, ~] = readmovielens('ml-data/', path);
        
        user_sex = user{1,3};
        
        females = getGenderIndex('female', user);
        males = getGenderIndex('male', user);
        
        young = getAgeIndex('young', user);
        adult = getAgeIndex('adult', user);
        elder = getAgeIndex('elder', user);
        
        % 1 = technician  6 = administrator 11 = entertainment 16 = marketing
        % 2 = programmer  7 = student       12 = librarian    17 = none
        % 3 = other       8 = lawyer        13 = homemaker    18 = healthcare
        % 4 = writer      9 = educator      14 = artist       19 = retired
        % 5 = executive  10 = scientist     15 = engineer
        educator_indx = getOccupationIndex('educator', user);
        programmer_indx = getOccupationIndex('programmer', user);
        technician_indx = getOccupationIndex('technician', user);
        writer_indx = getOccupationIndex('writer', user);
        executive_indx = getOccupationIndex('executive', user);
        administrator_indx = getOccupationIndex('administrator', user);
        student_indx = getOccupationIndex('student', user);
        lawyer_indx = getOccupationIndex('lawyer', user);
        scientist_indx = getOccupationIndex('scientist', user);
        entertainment_indx = getOccupationIndex('entertainment', user);
        librarian_indx = getOccupationIndex('librarian', user);
        homemaker_indx = getOccupationIndex('homemaker', user);
        artist_indx = getOccupationIndex('artist', user);
        engineer_indx = getOccupationIndex('engineer', user);
        marketing_indx = getOccupationIndex('marketing', user);
        none_indx = getOccupationIndex('none', user);
        healthcare_indx = getOccupationIndex('healthcare', user);
        retired_indx = getOccupationIndex('retired', user);
        
        % kmeans clustering, seeds to make them the same.
        rand( 'seed', 0 ); %#ok
        randn( 'seed', 0 ); %#ok
        
        % Makes a 2d plot based on the first two principal components. Running
        % kmeans clustering on it using k = 50 (default)
        
        if(length(k) == 1)
            if(strcmp(space, 'U'))
                cl = kmeans( U, k );
            elseif(strcmp(space, 'V'))
                cl = kmeans( V, k );
            end
        elseif(length(k) == 2)
           if(strcmp(space, 'U'))
                [cl k] = autofind( U, k );
            elseif(strcmp(space, 'V'))
                [cl k] = autofind( V, k );            
           end
        end
        
        % Create k colors using colorcube then convert to
        % rbg colors.
        cc = hsv2rgb(colorcube(k));
        
        % Get indices, remeber to change text accordingly.
        indx = females & (adult | young);
        
        % + 0.015 is just to position the text a bit above the points just plotted,
        % because of the additional LineWidth parameter added to scatter.
        figure(461)
        hold on
        if(strcmp(space, 'U'))
            for i=1:k
                plot( U_new(cl == i & indx,1), U_new(cl == i & indx,2), ...
                    '.', 'MarkerSize', 6, 'Color', cc(i, :) )
            end
        elseif(strcmp(space, 'V'))
            for i=1:k
                plot( V(cl == i,1), V(cl == i,2), '.', 'MarkerSize', 6, 'Color', cc(i, :) )
            end
        end
        if(strcmp(textOn, 'true'))
            text( U_new(females & young,1), U_new(females & young,2) + 0.001, strcat(user_sex(females & young),'y'), 'FontSize', 6);
            text( U_new(females & adult,1), U_new(females & adult,2) + 0.001, strcat(user_sex(females & adult),'a'), 'FontSize', 6);
        end
        
        xlabel( 'pc1' )
        ylabel( 'pc2' )
        title([ 'Showing space ' space ' using k = ' num2str(k) ' for coloring']);
        hold off
        
        indx = males & (young | adult);
        
        figure(462)
        hold on
        if(strcmp(space, 'U'))
            for i=1:k
                plot( U_new(cl == i & indx,1), U_new(cl == i & indx,2), ...
                    '.', 'MarkerSize', 6, 'Color', cc(i, :) )
            end
        elseif(strcmp(space, 'V'))
            for i=1:k
                plot( V(cl == i,1), V(cl == i,2), '.', 'MarkerSize', 6, 'Color', cc(i, :) )
            end
        end
        if(strcmp(textOn, 'true'))
            text( U_new(males & young,1), U_new(males & young,2) + 0.001, strcat(user_sex(males & young),'y'), 'FontSize', 6);
            text( U_new(males & adult,1), U_new(males & adult,2) + 0.001, strcat(user_sex(males & adult),'a'), 'FontSize', 6);
        end
        
        xlabel( 'pc1' )
        ylabel( 'pc2' )
        title([ 'Showing space ' space ' using k = ' num2str(k) ' for coloring']);
        hold off        
        
    end

    % Larger k's seems to get lower mean(sumd) atleast when dealing with
    % the latent user space.
    function [optimal_idx optimal_k] = autofind(X, kv)
        % Just find a k.
        optimal_k = 0;
        % Select high then optimize down.
        optimal_meansumd = 10e9;
        
        fprintf('Auto finding k in range: %d--%d\n', kv(1), kv(2))
        
        for k = kv(1):kv(2)
            [idx, ~, sumd ] = kmeans( X, k );
            if(mean(sumd) < optimal_meansumd)
                optimal_meansumd = mean(sumd);
                optimal_k = k;
                optimal_idx = idx;
                fprintf('Better k found, new k: %d, new sumd: %f\n', optimal_k, optimal_meansumd)
            end
        end
    end
                
    function run()
        addpath('./hfunc/')

        space = 'U';
%         space = 'V';
        k = [3];
        textOn = 'true';
        
        q46knn(k, space, textOn);
    end
    run()
end