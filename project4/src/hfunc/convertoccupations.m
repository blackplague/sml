function occupationmatrix = convertoccupations( occucell )
%convertoccupations Function for converting occupations to numerical
% values, this too uses uint8 to convert. This is standard ASCII
% convertion. We perform scaling of the resulting occupation dividing it by
% 10e35.
        
    occupationmatrix = zeros(size(occucell));
    
    for i = 1:size(occucell,1)
        tmp = uint8(occucell{i});
        string = '';
        for j=1:size(tmp,2)
            string = strcat(string, num2str(tmp(1,j)));
        end
        occupationmatrix(i,1) = str2num(string)/10e35;
    end
end