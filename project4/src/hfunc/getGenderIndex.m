function [ gender_indx ] = getGenderIndex( gender_string, user, path )
%getGenderIndex return indices based on user gender (Automatically lowers
%input string).
%   Accepts following strings:
%   male: return index of all males
%   female: return index of all females

    if(nargin > 2)
        path = 'u.data';
        [~, ~, user, ~, ~] = readmovielens('ml-data/', path );
    end
    
    user_sex = user{1,3};
    user_sex = convertsex(user_sex);

    switch lower(gender_string)
        case {'male'}
            gender_indx = user_sex == 0;
        case {'female'}
            gender_indx = user_sex == 1;
    end

end

