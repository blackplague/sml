function [ output_args ] = q2_2( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function q2_2main(Phi1, Phiu2, t)
       
        beta = 1;
        alpha = 2;
        
        % This is for eq 3.49 p(w|t)
        inv_SN1 = (1/alpha)*eye(length(w)) + beta*Phi1'*Phi1;  % eq. 3.53
        mN1 = beta*inv(inv_SN2)*Phi1'*t;                       % eq. 3.54

        inv_SN2 = (1/alpha)*eye(length(w)) + beta*Phi2'*Phi2;  % eq. 3.53
        mN2 = beta*inv(inv_SN2)*Phi2'*t;                       % eq. 3.54
        
        y1 = mN1*t
        
        y2 = mN2*t;
        
        
    end


    function run( ~ )

        selection1 = [ 4 7 8 9];
        
        selection2 = [ 8 ];
        
        [training_dataset test_dataset] = readbodyfat;
        
        [ Phi1 Phi2 ] = q2_1main(selection1, selection2, training_dataset, test_dataset);
        
        q2_2main(Phi1, Phi2, training_data)
        
    end

    run()

end

