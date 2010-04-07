function datematrix = convertdate( datecell )
%convertdate Function for converting release dates into numerical values if the
% release date is unknown it is mapped to 0 (zero)

    datematrix = zeros(size(datecell));

    for i = 1:size(datecell,1)
        if strcmp('', datecell{i})
            datematrix(i, 1) = 0;
        end
        datematrix(i, 1) = datenum(datecell{i});
    end
end