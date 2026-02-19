function AUC = ea_logit_regression(Ihat_train, Ihat, Improvement, training, test)

% Fit logit model, compute ROC and find the optimal threshold.
% Compute confustion matrix for the test set (can be the same as training)

% if Ihat_train was not provided, then we have in-sample analysis
if isstruct(Ihat_train)
    % imported model 
    mdl = Ihat_train.mdl;
    scores_thresh = Ihat_train.scores_thresh;
    AUC = nan;
else
    if Ihat_train == 0
        Ihat_train = Ihat(training);
    end

    % first, we fit a logit function for our binary prediction
    mdl = fitglm(Ihat_train,Improvement(training),'Distribution','binomial','Link','logit');
    
    % Ihat_train_fake = Ihat_train(Improvement(training) == 1,1);
    % % % FTG
    % % Ihat_train_balanced = [Ihat_train;Ihat_train_fake;Ihat_train_fake];
    % % Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
    % 
    % % % ERNA (Berlin 2mA and all)
    % Ihat_train_balanced = [Ihat_train;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
    % Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
    
    
    % ERNA (Cologne, 2mA)
    % Ihat_train_balanced = [Ihat_train;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
    % Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
    
    % % ERNA (Cologne, all)
    % Ihat_train_balanced = [Ihat_train;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
    % Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
    
    
    % % allTransientEntrain
    %Ihat_train_balanced = [Ihat_train;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
    %Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
    
    %Ihat_train_balanced = Ihat_train;
    %Improvement_balanced = Improvement(training);
    
    %mdl = fitglm(Ihat_train_balanced,Improvement_balanced ,'Distribution','binomial','Link','logit');
    
    % second, we run ROC curve analysis
    scores = mdl.Fitted.Probability;
    [X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);
    %[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement_balanced,scores,1);
    figure, plot(X,Y, 'k', 'linew', 1.5)
    set(gcf,'color','w');
    hold on
    plot(OPTROCPT(1),OPTROCPT(2),'ro', 'MarkerSize',10)
    xlabel('False positive rate') 
    ylabel('True positive rate')
    txt = ['AUC: ' num2str(AUC)];
    text(0.7,0.1,txt)
    title('ROC for Classification by Logistic Regression')
    
    % optimal threshold on the classifier
    scores_thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
    %scores_thresh = 0.5;
    
    % plot logit fit for training
    vec_val = min(Ihat_train):1:max(Ihat_train);
    figure
    set(gcf,'color','w');
    lims=[min(Ihat_train)-0.1*(max(Ihat_train)-min(Ihat_train)), max(Ihat_train)+0.1*(max(Ihat_train)-min(Ihat_train))];
    subplot(4,1,1)
    subtitle('Response');
    col=ea_color_wes('lifeaquatic');
    g=ea_raincloud_plot(Ihat_train(Improvement(training)==1)','color',col(3,:),'box_on',1);
    a1=gca;
    set(a1,'ytick',[])
    set(gca, 'xlim', lims)
    a1.YLabel.String='Response';
    a1.XLabel.String='Fiberscore';
    a1.Box='off';
    title('Logistic Regression for Training Cohort')
    
    subplot(4,1,[2 3])
    plot(vec_val', predict(mdl,vec_val'),'k', 'linew', 1.5)
    xlabel('Fiberscore'), ylabel('Response')
    set(gca, 'xlim', lims); box off
    %plot(vec_val', predict(mdl,vec_val'),Ihat_train,Improvement(training),'s')
    %plot(Ihat_av,Improvement,'s')
    
    subplot(4,1,4)
    subtitle('Control');
    col=ea_color_wes('lifeaquatic');
    g=ea_raincloud_plot(Ihat_train(Improvement(training)==0)','color',col(1,:),'box_on',1);
    a1=gca;
    set(a1,'ytick',[])
    set(gca, 'xlim', lims)
    a1.YLabel.String='Control';
    a1.XLabel.String='Fiberscore';
    a1.Box='off';
end


% prediction for test based on the logit model
scores_test = predict(mdl,Ihat(test));
Ihat_prediction = scores_test > scores_thresh;

% get the confussion matrix (this can be done on the test set now)
figure
cm = confusionchart(logical(Improvement(test)), Ihat_prediction);
set(gcf,'color','w');

tp = sum((Ihat_prediction == 1) & (Improvement(test) == 1));
fp = sum((Ihat_prediction == 1) & (Improvement(test) == 0));
tn = sum((Ihat_prediction == 0) & (Improvement(test) == 0));
fn = sum((Ihat_prediction == 0) & (Improvement(test) == 1));

sensitivity = tp/(tp + fn);  % TPR
specificity = tn/(tn + fp);  % TNR
precision = tp/(tp + fp);
f1 = 2*sensitivity*precision/(sensitivity+precision);

cm.Title = ['Sensitivity: ', sprintf('%.2f',sensitivity), '; ', 'Specificity: ', sprintf('%.2f',specificity), '; ', 'F1: ', sprintf('%.2f',f1)];


if (size(Improvement(test),1) == sum(Improvement(test))*2) || size(Improvement(test),1)==131
    % likely predicting a shell, do binomial test

    % Binomial test across the threshold
    % bad practice!
    
    % stupid way to define it
    n_trials = size(Improvement(test),1);
    k_successes = sum((Ihat_prediction & Improvement(test)) | (~Ihat_prediction & ~Improvement(test)));
    
    p_chance = 0.5;          % Probability of success by chance (0.5 for binary)
    alpha = 0.05;            % Significance level
    
    % 2. Calculate the p-value
    % We use 'upper' because we want to know if accuracy is GREATER than chance.
    % We subtract 1 from k because binocdf(k, n, p, 'upper') calculates P(X > k).
    % To get P(X >= k), we need P(X > k-1).
    p_value = binocdf(k_successes - 1, n_trials, p_chance, 'upper');
    
    % 3. Display results
    fprintf('--- Binomial Test Results ---\n');
    fprintf('Accuracy: %.2f%%\n', (k_successes/n_trials)*100);
    fprintf('p-value:  %.4f\n', p_value);
    
    if p_value < alpha
        fprintf('Result: Significant! (Reject Null Hypothesis)\n');
    else
        fprintf('Result: Not significant. (Fail to reject Null Hypothesis)\n');
    end
    
    % 4. Optional: Visualize the Null Distribution
    figure
    x = 0:n_trials;
    y = binopdf(x, n_trials, p_chance);
    bar(x, y, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
    hold on;
    stem(k_successes, binopdf(k_successes, n_trials, p_chance), 'r', 'LineWidth', 2);
    title('Binomial Distribution (Null Hypothesis)');
    xlabel('Number of Correct Predictions');
    ylabel('Probability');
    legend('Chance Distribution', 'Your Model');
    grid on;
end

[z_score, p_value] = delong_test_independent(Improvement(test), Ihat_prediction, Improvement(test), scores_test_amp);





end
