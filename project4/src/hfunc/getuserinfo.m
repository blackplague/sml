function [userfeatures ratingvector ratings] = getuserinfo( path )
%getuserinfo( path ) returns [ userfeatures ratings ]
    % X = user x movie rating
    
    [ ratings, timestamps, user, movie, genre ] = ...
        readmovielens('ml-data/', [ path ]);
    
    % movie_id 1..943
    movie_id = movie{1,1};
    % movie title : string
%     movie_title = movie{1,2};
    % movie release
%     movie_release = movie{1,3};
%     movie_videorelease = movie{1,4};
%     movie_url = movie{1,5};
    movie_genrecoded = movie{1,6};
%     movie_genrestring = movie{1,7};
    
    user_id = user{1,1};
    user_age = user{1,2};
    user_sex = user{1,3};
    user_occupation = user{1,4};
    user_zipcode = user{1,5};
    
    user_sex = convertsex( user_sex );
    user_zipcode = convertzipcodes( user_zipcode );
    user_occupation = convertoccupations( user_occupation );
    
    % x is preallocated for speed.
    userfeatures = zeros(length(find(ratings)), 26);
    
    ratingvector = zeros(length(find(ratings)),1);
    
    % Counter for making sure each user_id | movie_id entry has it own
    % row
    counter = 1;
    
    % Construction of x matrix
    % uid | mid | user_age | user_sex | user_zipcode | user_occupation
    % | movie_genre | time_stamp
    % This loop constructs the x_nm matrix, having
    % x_nm = [x_user,n x_item,m timestamp_nm]
    for n=1:size(ratings,1)
        for m=1:size(ratings,2)
            if ratings(n,m) ~= 0
                userfeatures(counter, :) = [ user_id(n) movie_id(m) ... 
                    user_age(n) user_sex(n) user_zipcode(n) ... 
                    user_occupation(n) movie_genrecoded(m,:) ...
                    timestamps(n,m)/10e8];
                ratingvector(counter,1) = ratings(n,m);
                counter = counter + 1;
            end
        end
    end
end