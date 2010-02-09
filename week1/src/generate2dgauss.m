function y = generate2dgauss(number)

    x = randn(2, number);

    mu = [1.0; 1.0];
    sigma = [0.3 0.2; 0.2 0.2];
    
    L = chol(sigma, 'lower');
    
    y = 1.0 + L * x;
        
end
