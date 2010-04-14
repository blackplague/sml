function p = q42( Theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    function [m1 m2 C1 C2 prior1 train_err reject_rate p] = q42trainclass( uf, dr, Theta )
        
        if nargin < 2
            Theta = 0;
        end

        class1_indx = dr == 1;
        class2_indx = dr == 2;
        
        N1 = sum(class1_indx);
        N2 = sum(class2_indx);
        N = N1+N2;
        
%         size(class1_indx)
%         size(uf)
        
        m1 = mean(uf(class1_indx, 3:end));
        m2 = mean(uf(class2_indx, 3:end));
        
        C1 = cov(uf(class1_indx, 3:end));
        C2 = cov(uf(class2_indx, 3:end));
        
        prior1 = N1/N;
        prior2 = 1 - prior1;
        
        pclass = zeros(length(dr), 2);
        p = zeros(length(dr), 2);
        
        for i = 1:length(dr)
            pclass(i,1) = mvnpdf( uf(i,3:end), m1, C1 )*prior1;
            pclass(i,2) = mvnpdf( uf(i,3:end), m2, C2 )*prior2;
        end
        
        p(:,1) = logsig(log10(pclass(:,1)./pclass(:,2)));
        p(:,2) = 1-p(:,1);
        
        pdr = predicteddiscreterating(p, Theta);
        
        [train_err reject_rate] = misclassificationrate(dr, pdr);
        
    end

    function [err_rate rej_rate] = misclassificationrate(dr, pdr)

        misclass = 0;        
        correctionNumber = 0;
        % Loop to find misclassified
        for i=1:length(dr)
            if(pdr(i) == -1)
                correctionNumber = correctionNumber + 1;
                continue;
            elseif(dr(i) ~= pdr(i))
                misclass = misclass + 1;
            end
        end
        
        err_rate = (misclass / (length(dr)-correctionNumber))*100;
        rej_rate = (correctionNumber / length(dr))*100;
    end

    function pdr = predicteddiscreterating(p, Theta)
        % Init all to something that not a class
        pdr = -1*ones(length(p),1);
        
        % Assign class by using the max of probability of the two classes
        for i=1:length(p)
            [prob, I] = max(p(i,:));
            if(prob < Theta)
                pdr(i) = -1;
            else
                pdr(i) = I;
            end
        end
    end

    function [test_err reject_rate] = classifier(m1, m2, C1, C2, prior, uf, dr, Theta)
        
        if nargin < 8
            Theta = 0;
        end
                
        prior1 = prior;
        prior2 = 1 - prior;
        
        pclass = zeros(length(dr), 2);
        p = zeros(length(dr), 2);
        
        for i = 1:length(dr)
                pclass(i,1) = mvnpdf( uf(i,3:end), m1, C1 )*prior1;
                pclass(i,2) = mvnpdf( uf(i,3:end), m2, C2 )*prior2;
        end

        p(:,1) = logsig(log10(pclass(:,1)./pclass(:,2)));        
        p(:,2) = 1 - p(:,1);
        
        pdr = predicteddiscreterating(p, Theta);
        
        [test_err reject_rate] = misclassificationrate(dr, pdr);
        
    end

    function [newUF newDR] = discretize( uf, rv )
    %discretize - Function for discretizing the rating matrix, 
    %it maps 1 and 2 -> 2 for 'dislike', 4 and 5 -> 1 for 'like' and 
    %just throws the rest out.        
        N = length(find(rv == 1 | rv == 2 | rv == 4 | rv == 5));
        
        newUF = zeros(N, size(uf,2));
        newDR = zeros(N,1);
        
        c = 1;
        for i=1:length(rv)
            switch rv(i,1)
                case {1,2}
                    newDR(c,1) = 2;
                    newUF(c,:) = uf(i,:);
                    c = c + 1;
                case {4,5}
                    newDR(c,1) = 1;
                    newUF(c,:) = uf(i,:);
                    c = c + 1;
            end
        end   
    end

    function makeTableTex(tr_err, te_err, Theta)

        fprintf('\\begin{table}[!htbp]\n')
        fprintf('\\centering\n')
        fprintf('\\begin{tabular}{r r r}\n')
        fprintf('Set number & Misclassification training & Misclassification test \\\\\n')
        fprintf('\\hline\n')
        for i = 1:length(tr_err)
            fprintf('%d & %2.5f \\%% & %2.5f \\%% \\\\\n', i, tr_err(i), te_err(i))
        end
        fprintf('Mean & %2.5f \\%% & %2.5f \\%% \\\\\n', mean(tr_err), mean(te_err))
        fprintf('\\end{tabular}\n')
        fprintf('\\caption{Shows misclassification rates for training and test, along with the mean with $\\Theta = %1.2f$.}\n', Theta)
        fprintf('\\label{tab:q42misctrtet%d}\n', Theta*100);
        fprintf('\\end{table}\n')        
    end

if nargin < 1
    Theta = 0;
end

    function p = run(Theta)
        
        addpath('./hfunc/')
        
        % Preload for speed.
        uf_train = cell(1,5); rv_train = cell(1,5); dr_train = cell(1,5);
        uf_test  = cell(1,5); rv_test  = cell(1,5); dr_test  = cell(1,5);
        for i = 1:length(uf_train)
            [ uf_train{i} rv_train{i}] = getuserinfo( ['u' num2str(i) '.base'] );
            [ uf_train{i} dr_train{i}] = discretize( uf_train{i}, rv_train{i});
            [ uf_test{i}  rv_test{i}]  = getuserinfo( ['u' num2str(i) '.test'] );
            [ uf_test{i}  dr_test{i}]  = discretize( uf_test{i}, rv_test{i});
        end
        
%         Theta = 0.5:0.05:1;

        % Init all the variables need to make this run.
        m1 = zeros(5,24); m2 = zeros(5,24); C1 = cell(5,1); C2 = cell(5,1);
        tr_err = zeros(5,1); prior = zeros(5,1); te_err = zeros(5,1);
        reject_train = zeros(5,1); reject_test = zeros(5,1);
        rr_tr = zeros(length(Theta)); rr_te = zeros(length(Theta));
        tre = zeros(length(Theta)); tee = zeros(length(Theta));        

        % Run a loop over Theta and the 5 training / test sets and
        % calculate the cross validation to be plottet later, along with
        % rejection rate.
        for j = 1:length(Theta)
            for i=1:size(tr_err,1)
                [m1(i,:) m2(i,:) C1{i} C2{i} prior(i) tr_err(i) reject_train(i) p] = ...
                    q42trainclass(uf_train{i}, dr_train{i}, Theta(j));
                fprintf('Done training %d, Theta %f\n', i, Theta(j))
                [te_err(i) reject_test(i)] = classifier(m1(i,:), ...
                    m2(i,:), C1{i}, C2{i}, prior(i), uf_test{i}, ...
                    dr_test{i}, Theta(j));
                fprintf('Done classifying %d, Theta %f\n', i, Theta(j))
            end
            tre(j) = tre(j) + mean(tr_err);
            tee(j) = tee(j) + mean(te_err);
            rr_tr(j) = rr_tr(j) + mean(reject_train);
            rr_te(j) = rr_te(j) + mean(reject_test);
        end

        tre(:) = tre(:)/length(tr_err);
        tee(:) = tee(:)/length(te_err);
        rr_tr(:) = rr_tr(:)/length(tr_err);
        rr_te(:) = rr_te(:)/length(te_err);
        
        fprintf('\n')
        
        makeTableTex(tr_err, te_err, Theta);            
        
%         fig1 = figure(100);
%         hold on
%         plot(Theta, tre, 'b--')
%         plot(Theta, tee, 'r-')
%         xlim([min(Theta); max(Theta)])
%         ylim([min(min(tre,tee))-10; max(max(tre,tee))+10])
%         xlabel('Threshold (\Theta)')
%         ylabel('Misclassification (%)')
%         title('Training and test misclassification as a function of Theta')
%         legend('Training', 'Test')
%         hold off
%         saveas(fig1, '../report/images/q42miscft', 'epsc')
% 
%         fig2 = figure(110);
%         hold on
%         plot(Theta, tre, 'b--')
%         plot(Theta, tee, 'r-')
%         xlim([min(Theta); max(Theta)])
%         ylim([min(min(rr_tr,rr_te))-10; max(max(rr_tr,rr_te))+10])
%         xlabel('\Theta')
%         ylabel('Rejction rate (%)')
%         title('Training and test rejection rate as a function of Theta')
%         legend('Training', 'Test')
%         hold off
%         saveas(fig2, '../report/images/q42rrft', 'epsc')        
        
    end
    p = run(Theta);
end