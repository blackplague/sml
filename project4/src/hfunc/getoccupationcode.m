function [ ret_str ] = getoccupationcode( str )
%getoccupationcdoe Takes string and returns occupation code
%   Detailed explanation goes here

        tmp = uint8(str);
        string = '';
        for j=1:size(tmp,2)
            string = strcat(string, num2str(tmp(1,j)));
        end
        ret_str = str2num(string)/10e35;
end

