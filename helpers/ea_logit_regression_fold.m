function [Ihat_prediction,scores_test,log_loss,AUC] = ea_logit_regression_fold(Ihat_train, Ihat, Improvement, training, test)

% Fit logit model, compute ROC and find the optimal threshold.
% Compute confustion matrix for the test set (can be the same as training)

% if Ihat_train was not provided, then we have in-sample analysis
if Ihat_train == 0
    Ihat_train = Ihat(training);
end

% first, we fit a logit function for our binary prediction
mdl = fitglm(Ihat_train,Improvement(training),'Distribution','binomial','Link','logit');


% let's upsample response outcomes by copying the corresponing Ihat_train
% and Improvement entries
% % ERNA (Berlin all and 2mA)
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];


% % ERNA (Cologne 2mA)
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];

% % ERNA (Berlin and Cologne, 2mA)
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];



% mdl = fitglm(Ihat_train_balanced,Improvement_balanced ,'Distribution','binomial','Link','logit');
% 
% % FTG
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];

% allTransientEntrainAdjusted
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];

% allTransientEntrain
% Ihat_train_fake = Ihat_train(1,Improvement(training) == 1)';
% Ihat_train_balanced = [Ihat_train';Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake;Ihat_train_fake];
% Improvement_balanced = [Improvement(training);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1);true(size(Ihat_train_fake,1),1)];
% 


%mdl = fitglm(Ihat_train_balanced,Improvement_balanced ,'Distribution','binomial','Link','logit');



% second, we run ROC curve analysis
scores = mdl.Fitted.Probability;
%[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement_balanced,scores,1);
[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);

% optimal threshold on the classifier
scores_thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
%scores_thresh = 0.5;

% prediction for test based on the logit model
scores_test = predict(mdl,Ihat(test));
Ihat_prediction = scores_test > scores_thresh;

% Compute Log Loss
log_loss = -mean(Improvement(test) .* log(scores_test) + (1 - Improvement(test)) .* log(1 - scores_test));


end
