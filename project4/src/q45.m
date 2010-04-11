function [ output_args ] = q45( Comment )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function q45pca_a( X, Comment )
        addpath('hfunc/')
        
        path = 'u.data';
        
        Y = X - repmat( mean( X ), size( X, 1 ), 1 );
        
        % running SVD/PCA
        [ U S V ] = svd( Y, 0 );
        
        [ratings, ~, user, ~, ~] = readmovielens('ml-data/', path);
        
        user_age = user{1,2};
        user_sex = user{1,3};
        user_occupation = user{1,4};
        user_zipcode = user{1,5};
        
        r_m50_dis = ratings(:,50) < 3 & ratings(:,50) ~= 0;
        r_m50_like = ratings(:,50) > 3;
        
        uz = convertzipcodes( user_zipcode );
        uo = convertoccupations( user_occupation );

        user_sex = convertsex(user_sex);
        
        % 1 = technician  6 = administrator 11 = enterainment 16 = marketing      
        % 2 = programmer  7 = student       12 = librarian    17 = none
        % 3 = other       8 = lawyer        13 = homemaker    18 = healthcare
        % 4 = writer      9 = educator      14 = artist       19 = retired
        % 5 = executive  10 = scientist     15 = engineer
        useroccupations = {'technician'; 'programmer'; 'other'; 'writer'; ...
            'executive'; 'administrator'; 'student'; 'lawyer'; 'educator'; ...
            'scientist'; 'enterainment'; 'librarian'; 'homemaker'; ...
            'artist'; 'engineer'; 'marketing'; 'none'; 'healthcare'; 'retired'};
        
        selected_uo1 = 19;
        selected_uo2 = 18;
        selected_uo3 = 13;
        
        % Education
        uo_idx1 = uo == getoccupationcode(useroccupations{selected_uo1});
        uo_idx2 = uo == getoccupationcode(useroccupations{selected_uo2});
        uo_idx3 = uo == getoccupationcode(useroccupations{selected_uo3});
        
        % Sexes
        females = user_sex == 1;
        males = user_sex == 0;
        
        % Ages
        young = user_age < 24;
        adult = 25 < user_age & user_age < 49;
        elder = 50 < user_age;
        
        if(strcmp(Comment,'true'))
            fprintf('\nSize of U: %d x %d\n', size(U,1), size(U,2))
            fprintf('Size of S: %d x %d\n', size(S,1), size(S,2))
            fprintf('Size of V: %d x %d\n', size(V,1), size(V,2))
            fprintf('Size of User: %d x %d\n', size(User,1),size(User,2))
            fprintf('Size of Movie: %d x %d\n', size(Movie,1), size(Movie,2))        
        end
        
        % Plotting users based on sexes, female is red and male is blue
        figure(1)
        hold on
        plot(U(females,1), U(females,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
        plot(U(males,1), U(males,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3)
        legend('Female', 'Male')
        title('PCs - user sexes')
        hold off

%         % Plotting users based on their taste towards movie_id 50
%         figure(2)
%         hold on
%         plot(U(r_m50_like,1), U(r_m50_like,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(r_m50_dis,1), U(r_m50_dis,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend('Like', 'Dislike')
%         title('PCs Like / Dislike')
%         hold off
% 
%         % Plotting users based on occupations
%         figure(3)
%         hold on
%         plot(U(uo_idx1,1), U(uo_idx1,2), 'y.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(uo_idx2,1), U(uo_idx2,2), 'm.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(uo_idx3,1), U(uo_idx3,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend(useroccupations{selected_uo1}, useroccupations{selected_uo2}, useroccupations{selected_uo3})
%         title('PCs occupations')
%         hold off
%         
%         % Plotting users based on their age.
%         figure(4)
%         hold on
%         plot(U(young,1), U(young,2), 'y.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(adult,1), U(adult,2), 'm.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(elder,1), U(elder,2), 'k.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend('Young age < 24','Adult 25 < age < 49','Elder 50 < age')
%         title('PCs age')
%         hold off
        
    end

    function q45pca_b( X, Comment )
        addpath('hfunc/')
        
        path = 'u.data';
        
        Y = X - repmat( mean( X ), size( X, 1 ), 1 );
        
        % running SVD/PCA
        [ U S V ] = svd( Y, 0 );
        
        [ratings, ~, user, ~, ~] = readmovielens('ml-data/', path);
        
        user_age = user{1,2};
        user_sex = user{1,3};
        user_occupation = user{1,4};
        user_zipcode = user{1,5};
        
        r_m50_dis = ratings(:,50) < 3 & ratings(:,50) ~= 0;
        r_m50_like = ratings(:,50) > 3;
        
        uz = convertzipcodes( user_zipcode );
        uo = convertoccupations( user_occupation );
        
        user_sex = convertsex(user_sex);
        
        % 1 = technician  6 = administrator 11 = enterainment 16 = marketing      
        % 2 = programmer  7 = student       12 = librarian    17 = none
        % 3 = other       8 = lawyer        13 = homemaker    18 = healthcare
        % 4 = writer      9 = educator      14 = artist       19 = retired
        % 5 = executive  10 = scientist     15 = engineer
        useroccupations = {'technician'; 'programmer'; 'other'; 'writer'; ...
            'executive'; 'administrator'; 'student'; 'lawyer'; 'educator'; ...
            'scientist'; 'enterainment'; 'librarian'; 'homemaker'; ...
            'artist'; 'engineer'; 'marketing'; 'none'; 'healthcare'; 'retired'};
        
        selected_uo1 = 19;
        selected_uo2 = 18;
        selected_uo3 = 13;
        
        % Education
        uo_idx1 = uo == getoccupationcode(useroccupations{selected_uo1});
        uo_idx2 = uo == getoccupationcode(useroccupations{selected_uo2});
        uo_idx3 = uo == getoccupationcode(useroccupations{selected_uo3});
        
        % Sexes
        females = user_sex == 1;
        males = user_sex == 0;
        
        % Ages
        young = user_age < 24;
        adult = 25 < user_age & user_age < 49;
        elder = 50 < user_age;
        
        if(strcmp(Comment,'true'))
            fprintf('\nSize of U: %d x %d\n', size(U,1), size(U,2))
            fprintf('Size of S: %d x %d\n', size(S,1), size(S,2))
            fprintf('Size of V: %d x %d\n', size(V,1), size(V,2))
            fprintf('Size of User: %d x %d\n', size(User,1),size(User,2))
            fprintf('Size of Movie: %d x %d\n', size(Movie,1), size(Movie,2))        
        end
        
        % Plotting users according to sexes.
%         fig1 = figure(1);
%         hold on
%         plot( X(females,:)*V(:,1), X(females,:)*V(:,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3  )
%         plot( X(males,:)*V(:,1), X(males,:)*V(:,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3  )
%         legend('Female', 'Male')
%         title('Based X*V')
%         hold off
% %         hold on
% %         text( V(:,1), V(:,2) + 0.01, cancer_factor, 'FontSize', 8 )
% %         hold off
%         xlabel( 'pc1' )
%         ylabel( 'pc2' )

        % Plotting users based on sexes, female is red and male is blue
        figure(11)
        hold on
        plot(U(females,1), U(females,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
        plot(U(males,1), U(males,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3)
        legend('Female', 'Male')
        title('PCs user sex')
        hold off

%         figure(12)
%         hold on
%         plot(U(r_m50_like,1), U(r_m50_like,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(r_m50_dis,1), U(r_m50_dis,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend('Like', 'Dislike')
%         title('PCs Like / Dislike (Movie id 50)')
%         hold off
% 
%         figure(13)
%         hold on
%         plot(U(uo_idx1,1), U(uo_idx1,2), 'y.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(uo_idx2,1), U(uo_idx2,2), 'm.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(uo_idx3,1), U(uo_idx3,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend(useroccupations{selected_uo1}, useroccupations{selected_uo2}, useroccupations{selected_uo3})
%         title('PCs occupations')
%         hold off
%         
%         figure(14)
%         hold on
%         plot(U(young,1), U(young,2), 'y.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(adult,1), U(adult,2), 'm.', 'MarkerSize', 6, 'LineWidth', 3)
%         plot(U(elder,1), U(elder,2), 'k.', 'MarkerSize', 6, 'LineWidth', 3)
%         legend('Young age < 24','Adult 25 < age < 49','Elder 50 < age')
%         title('PCs age')
%         hold off
        
    end

    function q45pca_c( X )
        addpath('hfunc/')
        
        path = 'u.data';
        texton = 'true';
        
        Y = X - repmat( mean( X ), size( X, 1 ), 1 );
        
        % running SVD/PCA
        [ U S V ] = svd( Y, 0 );
        
        [ratings, timestamps, user, movie, genre] = readmovielens('ml-data/', path);
        
        movie_id = movie{1,1};
        movie_title = movie{1,2};
        movie_genre = movie{1,7};
        
        % Mean vector 1 x 1682, contains mean rating for each movie.
        m_vec = meannonzero( ratings );
        % Rating vector 1 x 1682, contains how many have rated a movie.
        r_vec = numberofratings( ratings );
        
        bad   = m_vec < 1.5;
        poor  = 1.5 < m_vec & m_vec < 2.5;
        mid   = 2.5 < m_vec & m_vec < 3.5;
        good  = 3.5 < m_vec & m_vec < 4.5;
        super = 4.5 < m_vec;
        fprintf('\n')
        fprintf('Number of bad movies: %d\n', length(find(bad)))
        fprintf('Number of poor movies: %d\n', length(find(poor)))
        fprintf('Number of mid movies: %d\n', length(find(mid)))
        fprintf('Number of good movies: %d\n', length(find(good)))
        fprintf('Number of super movies: %d\n', length(find(super)))        
        
        % Plotting users according to sexes.

        fig23 = figure(23);
        hold on
        plot( U(bad,1), U(bad,2), 'k.', 'MarkerSize', 6, 'LineWidth', 3  )
        plot( U(poor,1), U(poor,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3  )
        plot( U(mid,1), U(mid,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3  )
        plot( U(good,1), U(good,2), 'm.', 'MarkerSize', 6, 'LineWidth', 3  )
        plot( U(super,1), U(super,2), 'c.', 'MarkerSize', 6, 'LineWidth', 3  )
        title('U')
        legend('Bad - Mean < 1.5', 'Poor - 1.5 < Mean < 2.5', ...
            'Mid - 2.5 < Mean < 3.5', 'Good - 3.5 < Mean < 4.5', ...
            'Super - 4.5 < Mean')
        xlabel( 'pc1' )
        ylabel( 'pc2' )        
        hold off
        saveas(fig23, ['../report/images/q45_' num2str(23)], 'epsc')
        
        fig24 = figure(24);
        hold on
        plot( U(super,1), U(super,2), 'r.', 'MarkerSize', 6, 'LineWidth', 3  )
        title('Top rated - Movies with a rating > 4.5')
        if(strcmp(texton, 'true'))
            text( U(super,1), U(super,2) + 0.001, movie_title(super), 'FontSize', 6 )
            text( U(super,1), U(super,2) + 0.002, movie_genre(super), 'FontSize', 6 )
            text( U(super,1), U(super,2) + 0.003, r_vec(super), 'FontSize', 6 )
        end
        xlabel( 'pc1' )
        ylabel( 'pc2' )        
        hold off
        saveas(fig24, ['../report/images/q45_' num2str(24)], 'epsc')

        starwars = 50;
        nearest_indx = findnearest(10, starwars, U(:,1), U(:,2));

        figure(25)
        hold on
        plot( U(:,1), U(:,2), 'b.', 'MarkerSize', 6, 'LineWidth', 3  )
        plot( U(nearest_indx,1), U(nearest_indx,2), 'mo', 'MarkerSize', 4, 'MarkerFaceColor', 'm'  )
        plot( U(starwars,1), U(starwars,2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'  )
        title('10 Movies nearest to Starwars (1977)')
        if(strcmp(texton, 'true'))
            text( U(nearest_indx,1), U(nearest_indx,2) + 0.001, movie_title(nearest_indx), 'FontSize', 8 )
            text( U(nearest_indx,1), U(nearest_indx,2) + 0.0016, movie_genre(nearest_indx), 'FontSize', 6 )
            text( U(starwars,1), U(starwars,2) + 0.001, movie_title(starwars), 'FontSize', 8 )
            text( U(starwars,1), U(starwars,2) + 0.0016, movie_genre(starwars), 'FontSize', 6 )
        end
        xlabel( 'pc1' )
        ylabel( 'pc2' )        
        hold off        
        
    end

    function run(Comment)
        
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
       
        if(strcmp(Comment, 'true'))
            fprintf('\n')
            fprintf('Size of Y: %d x %d\n', size(Y,1), size(Y,2))
            fprintf('Size of U: %d x %d\n', size(U,1), size(U,2))
            fprintf('Size of V: %d x %d\n', size(V,1), size(V,2))
        end
        
%         q45pca_a(Y_mf, Comment)
        
%         q45pca_b(U, Comment)
        
        q45pca_c( V )
        
    end

run(Comment)

end