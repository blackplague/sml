function [ output_args ] = q2_1( input_args )

    function [ Phi1 Phi2 ] = q2_1main(extract1, extract2, training_data, test_data)

        size_data1 = size(training_data);
        t = training_data(:,2);
        
        extracted_values1 = zeros(size_data1(1), length(extract1));
        extracted_values2 = zeros(size_data1(1), length(extract2));
        
        size_test1 = size(test_data);
        tt = test_data(:,2);
        
        extract_test1 = zeros(size_test1(1), length(extract1));
        extract_test2 = zeros(size_test1(1), length(extract2));
        
        for i=1:length(extract1)
            extracted_values1(:,i) = training_data(:,extract1(i));
        end
        
        for i=1:length(extract2)
            extracted_values2(:,i) = training_data(:,extract2(i));
        end
        
        for i=1:length(extract1)
            extract_test1(:,i) = test_data(:,extract1(i));
        end
        
        for i=1:length(extract2)
            extract_test2(:,i) = test_data(:,extract2(i));
        end
        
        size_ev1 = size(extracted_values1);
        size_ev2 = size(extracted_values2);
        size_tv1 = size(extract_test1);
        size_tv2 = size(extract_test2);
        phi1 = zeros(size_ev1(2)+1, size_ev1(1));
        phi2 = zeros(size_ev2(2)+1, size_ev2(1));
        phi3 = zeros(size_tv1(2)+1, size_tv1(1));
        phi4 = zeros(size_tv2(2)+1, size_tv2(1));
        
        phi1(1, :) = 1;
        phi1(2:end, :) = extracted_values1';
        phi2(1, :) = 1;
        phi2(2:end, :) = extracted_values2';
        
        phi3(1, :) = 1;
        phi3(2:end, :) = extract_test1';
        phi4(1, :) = 1;
        phi4(2:end, :) = extract_test2';
        
        phi1 = phi1';
        phi2 = phi2';

        phi3 = phi3';
        phi4 = phi4';
        
        Phi1 = pinv(phi1);
        Phi2 = pinv(phi2);
%         fprintf('\rSize Theta1: %d\nSize Theta2: %d', size(Theta1), size(Theta2))
        
        w_ml1 = Phi1*t;
        w_ml2 = Phi2*t;
%         fprintf('\rSize w_ml2: %d x %d', size(w_ml2,1), size(w_ml2,2))
%         fprintf('\rSize tt: %d x %d', size(tt,1), size(tt,2))
        
        result1 = phi1*w_ml1;
        result2 = phi2*w_ml2;

        result3 = phi3*w_ml1;
        result4 = phi4*w_ml2;
        
        vector1 = t - result1;
        vector2 = t - result2;
        
        vector3 = tt - result3;
        vector4 = tt - result4;
        
        RMS1 = sqrt((vector1'*vector1)/length(vector1));
        RMS2 = sqrt((vector2'*vector2)/length(vector2));
        RMS3 = sqrt((vector3'*vector3)/length(vector3));
        RMS4 = sqrt((vector4'*vector4)/length(vector4));

        fprintf('\rTraining RMS1 = %d',RMS1);
        fprintf('\rTest RMS1 = \t%d',RMS3);
        fprintf('\rTraining RMS2 = %d',RMS2);
        fprintf('\rTest RMS2 = \t%d',RMS4);        
        fprintf('\rDifference betweeen Training and Test, selection1: %d', RMS1 - RMS3)
        fprintf('\rDifference betweeen Training and Test, selection2: %d\n', RMS2 - RMS4)        
        
    end

    function run( ~ )
       
        selection1 = [ 4 7 8 9];
        
        selection2 = [ 8 ];
        
        [training_dataset test_dataset] = readbodyfat;
        
        q2_1main(selection1, selection2, training_dataset, test_dataset);
        
    end

    run()

end