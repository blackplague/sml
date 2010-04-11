function [ indx_vector ] = getOccupationIndex( occupation_string, user, path )
%getOccupationIndex return index vector given occupation_string
    % Valid occupation strings are:
    % 1 = technician  6 = administrator 11 = enterainment 16 = marketing
    % 2 = programmer  7 = student       12 = librarian    17 = none
    % 3 = other       8 = lawyer        13 = homemaker    18 = healthcare
    % 4 = writer      9 = educator      14 = artist       19 = retired
    % 5 = executive  10 = scientist     15 = engineer

    if(nargin > 2)
        path = 'u.data';
        [~, ~, user, ~, ~] = readmovielens('ml-data/', path );
    end
    
    user_occupation = user{1,4};
    uo = convertoccupations( user_occupation );

    indx_vector = uo == getoccupationcode(occupation_string);

end

