function [ output_args ] = q410( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    function sqrtDiffVector = powerCalculation(d, T, init)
    %powerCalculation - Accept d - damping factor, T - Transition matrix, 
    % init - initialization of p. init can be the following 'uni', 
    % for p_i^(0) = 1/N, 'top' for p(maxindex(p_true)) = 1, 
    % rest p(rest) = 0, 'bottom' same as above, just the lowest and 'randn'
    % to get p initialized with pseudo random normally distributed values,
    % which afterwards are normalized to sum to 1.
        p_true = q49(d);
        N = size(T,2);
        
        if(strcmpi(init, 'uni'))
            % Init p_i^(0) = 1/N
            p = ones(N,1)/N;
        elseif(strcmpi(init, 'top'))
            % Init p_i = 1 for top rank and p_j = 0 for j != i
            p = zeros(N,1);
            [~, idx] = max(p_true);
            p(idx) = 1;
        elseif(strcmpi(init, 'bottom'))
            % Init p_i = 1 for buttom rank and p_j = 0 for j != i
            p = zeros(N,1);
            [~, idx] = min(p_true);
            p(idx) = 1;
        elseif(strcmpi(init, 'randn'))
            % Init p_i to random values, then scales to sum to one
            p = randn(N,1);
            p = p/sum(p);
        end        
        
        sqrtDiff = 10e10;
        sqrtDiffVector = [10^10 10^5];
        
        counter = 1;
        while(sqrtDiff > 10^(-10))
            p(:,counter+1) = (1-d)/N*ones(N,1) + d*T*p(:,counter);
            sqrtDiff = sum((p(:,counter+1) - p_true).^2);
            sqrtDiffVector(counter) = sqrtDiff;
            fprintf('Iteration %d sqrtDiff: %2.11f\n', counter, sqrtDiff)
            counter = counter + 1;
            if(abs(sqrtDiffVector(end)-sqrtDiffVector(end-1)) < 10^(-11))
                break;
            end
        end        

    end

    function plotSqrtDiff( sqrtDiffVector, name, d )
        
        fig = figure(1);
        plot(1:length(sqrtDiffVector), sqrtDiffVector(:), 'b-')
        title(['Squared difference versus iterations (' name ') d = ' num2str(d)])
        xlabel('Iterations')
        ylabel('Squared difference')
        tmp = num2str(d);
        saveas(fig, ['../report/images/q410' name 'd_' tmp(end-1) tmp(end)], 'epsc')
        
    end
        
    function run()
        
        addpath('./hfunc/')
        
        T = createTransitionMatrix();
        d = 0.85; init = 'uni';

        sqrtDiffVector = powerCalculation(d, T, init);
                
        plotSqrtDiff( sqrtDiffVector, init, d)
    end

run()
end