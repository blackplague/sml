function [ output_args ] = getrandomuser( no )
%getrandomuser return no of users
%   Detailed explanation goes here

    rand('seed',42)
    
    unidrnd(943, 1, no)


end

