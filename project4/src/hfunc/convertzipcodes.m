function zipmatrix = convertzipcodes( zipcell )
%convertzipcodes Function for converting zip codes into numerical values, it uses
% uint8 to do this. This is standard ASCII convertion.
% We perform scaling of the zip codes dividing it by 10e8
        
    zipmatrix = zeros(size(zipcell));
        
    for i = 1:size(zipcell,1)
        tmp = uint8(zipcell{i});
        string = '';
        for j=1:size(tmp,2)
            string = strcat(string, num2str(tmp(1,j)));
            zipmatrix(i,1) = str2num(string)/10e8;
        end
    end
end