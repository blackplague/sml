function [ age_indx ] = getAgeIndex( age_string, user, path )
%getAgeIndex returns index of user within age groups
%   Accepts following string (Automatically lowers input string):
%   young: user age < 24
%   adult: 25 < user age < 49
%   elder: 50 < user age

    if (nargin > 2)
        path = 'u.data';
        [~, ~, user, ~, ~] = readmovielens('ml-data/', path );
    end
    
    user_age = user{1,2};
    
    switch lower(age_string)
        case {'young'}
            age_indx = user_age < 24;
        case {'adult'}
            age_indx = 25 < user_age & user_age < 49;     
        case {'elder'}
            age_indx = 50 < user_age;
    end
end

