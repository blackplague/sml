function [ output_args ] = q2_4( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function q2_4main(x1, x2)
       
        
        
        N1 = size(x1, 1);
        N2 = size(x2, 1);
                
        mu1 = mean(x1', 2);
        mu2 = mean(x2', 2);
        
        S1 = cov(x1)*((N1-1)/N1)
        S2 = cov(x2)*((N2-1)/N2)
        
        S = S1+S2
        
        plotgrid(x1, x2)
        
%         S1 = 1/N1 * (x1 - repmat(mu1', length(x1), 1))*(x1 - repmat(mu1', length(x1), 1))';
%         S2 = 1/N2 * (x2 - repmat(mu2', length(x2), 1))*(x2 - repmat(mu2', length(x2), 1))';
%         
%         size(S1), size(S2)
%         
%         N = N1 + N2;
%         
%         S = N1/N*S1 + N2/N*S2
        
    end


    function run()
       
        N1=50;
        N2=250;
        x1 = getdataconddensity(N1,1);
        x2 = getdataconddensity(N2,2);
%         t = 
        
        q2_4main(x1, x2)
        
    end

    run()

end

