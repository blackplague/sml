function [p pageRank] = q49( d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function p_inf = direct(d, T)

        N = size(T,2);
        p_inf = ((1-d)/N)*(eye(N)-d*T)^(-1)*ones(N,1);
    end

    function [ nUsers ] = getUsers(p, number, location)
        
        if nargin < 2
            number = 10;
        end
        if nargin < 3
            location = 'top';
        end        
        
        userSize = 943;
        p_users = p(1:userSize);

        if(strcmpi(location, 'top'))
            [pageRank, idx] = sort(p_users, 1, 'descend');
        elseif(strcmpi(location, 'bottom'))
            [pageRank, idx] = sort(p_users, 1, 'descend');
        end
        
        if(strcmpi(location, 'top'))
            nUsers = [idx(1:number), pageRank(1:number)];
        elseif(strcmpi(location, 'bottom'))
            nUsers = [idx(end-number+1:end), pageRank(end-number+1:end)];
        end
            
    end

    function topMovies = getMovies(p, number, location)
        
        if nargin < 2
            number = 10;
        end        
        if nargin < 3
            location = 'top';
        end
        
        userSize = 943;
        p_movies = p(userSize+1:end);
        
        if(strcmpi(location, 'top'))
            [pageRank, idx] = sort(p_movies, 1, 'descend');
        elseif(strcmpi(location, 'bottom'))
            [pageRank, idx] = sort(p_movies, 1);
        end
        topMovies = [idx(1:number), pageRank(1:number)];
    end

%     function printPageRank( printThis ) 
%         
%     end
        

    function makeUserTable( users, d, U )
       
        Uid = U{1,1};
        Uage = U{1,2};
        Usex = U{1,3};
        Uoccupation = U{1,4};
        Uzipcode = U{1,5};
        
        idx = users(:,1);
        Uid = Uid(idx);
        Uage = Uage(idx);
        Usex = Usex(idx);
        Uoccupation = Uoccupation(idx);
        Uzipcode = Uzipcode(idx);
        
        fprintf('\\begin{table}[!htbp]\n')
        fprintf('\\centering\n')
        fprintf('\\begin{tabular}{llllll}\n')
        fprintf('User id & Age & Gender & Occupation & Zipcode & PageRank \\\\\n')
        fprintf('\\hline\n')
        for i=1:length(users)
            fprintf('%d & %d & %s & %s & %s & %2.6f \\\\\n', ...
                Uid(i), Uage(i), Usex{i}, ...
                Uoccupation{i}, Uzipcode{i}, users(i,2))
        end
        fprintf('\\end{tabular}\n') 
        fprintf('\\caption{X %d users according to PageRank, using $d = %1.2f$}\n', length(users), d)
        fprintf('\\label{tab:q49[X]ud%d}\n', d*100)
        fprintf('\\end{table}\n')
    end

    function makeMovieTable( movies, d, M, R )
        
        Mid = M{1,1};
        Mtitle = M{1,2};
        Mgenre = M{1,7};

        idx = movies(:,1);
        Mid = Mid(idx);
        Mtitle = Mtitle(idx);
        Mgenre = Mgenre(idx);
        meanMovieRating = meannonzero(R);        
        
        fprintf('\\begin{table}[!htbp]\n')
        fprintf('\\hspace{-2cm}\n')
        fprintf('\\begin{tabular}{lllll}\n')
        fprintf('Movie id & Title (Year) & Genre & Mean rating & PageRank \\\\\n')
        fprintf('\\hline\n')
        for i=1:length(movies)
            fprintf('%d & %s & %s & %2.3f & %2.6f \\\\\n', ...
                Mid(i), Mtitle{i}, Mgenre{i}, ...
                meanMovieRating(idx(i)), movies(i,2))
        end
        fprintf('\\end{tabular}\n') 
        fprintf('\\caption{X %d movies according to PageRank, using $d = %1.2f$}\n', length(movies), d)
        fprintf('\\label{tab:q49[X]md%d}\n', d*100)
        fprintf('\\end{table}\n')
        
    end

if nargin < 1
    d = 0.85;
end

    function [p pageRank] = run(d)
        
        addpath('./hfunc/')
        
        T = createTransitionMatrix();
        
%         n = 10;
        
        d=0.10:0.01:0.95;
        pageRank = zeros(size(d,2),2);
        n=1;
        
        for i = 1:length(d)
            p = direct(d(i), T);
            pageRank(i,:) = getUsers(p,n,'top');
        end
        
        fig = figure(2);
        plot(d, pageRank(:,2), 'b-')
        title('PageRank for top scoring user as a function of damping factor')
        xlabel('Damping factor')
        ylabel('PageRank')
        saveas(fig, '../report/images/q49topuserd', 'epsc')
        
%         usersTop = getUsers(p, n, 'top');
%         moviesTop = getMovies(p, n, 'top');
%         usersBottom = getUsers(p, n, 'bottom');
%         moviesBottom = getMovies(p, n, 'bottom');
       
%         [ R, ~, U, M, ~] = readmovielens('ml-data/', 'u.data');
%         
%         makeUserTable(usersTop, d, U);
%         makeUserTable(usersBottom, d, U);
%         makeMovieTable(moviesTop, d, M, R);
%         makeMovieTable(moviesBottom, d, M, R);
    end
[p pageRank] = run(d);
end