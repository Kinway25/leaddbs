function AUC = ea_logit_regression(Ihat_train, Ihat, Improvement, training, test, gpatsel, obj)

% Fit logit model, compute ROC and find the optimal threshold.
% Compute confustion matrix for the test set (can be the same as training)

disp_regression_plots = 1;

% if Ihat_train was not provided, then we have in-sample analysis
if Ihat_train == 0
    Ihat_train = Ihat(training);
end


% % % let's iterate only over training and create a table
% I_training = Improvement(training);
% for pt_j = 1:length(gpatsel)
% 
%     if training(pt_j) == 1
%         temp2=strsplit(obj.M.patient.list{gpatsel(pt_j)},'_');
%         pt_ID2 = temp2{2}; 
%         %trajDepth2 = temp2{end}(4:end);
%         side2 = temp2{3};
%         side_label = side2;
%         
%         if contains(obj.M.patient.list{gpatsel(pt_j)}, '_fl.')
%             amp = str2num(temp2{end-1}(1:3));
%             trajDepth2 = temp2{end-1}(4:end);
%         else
%             amp = str2num(temp2{end}(1:3));
%             trajDepth2 = temp2{end}(4:end);
%         end
% 
%         row_to_add = {pt_ID2,trajDepth2,side_label,Ihat(pt_j),amp, Improvement(pt_j)};
%         if ~exist('PatientsIntraOp_FF','var')
%             PatientsIntraOp_FF = table;
%             PatientsIntraOp_FF.Code = pt_ID2;
%             PatientsIntraOp_FF.TrajDepth = trajDepth2;
%             PatientsIntraOp_FF.Side = side_label;
%             PatientsIntraOp_FF.FF_score = Ihat(pt_j);
%             PatientsIntraOp_FF.ClinAny = amp;
%             PatientsIntraOp_FF.Response = Improvement(pt_j);
%         else
%             PatientsIntraOp_FF = [PatientsIntraOp_FF;row_to_add];
%         end
%     end
% end
% 
% PatientsIntraOp_FF.Code = categorical(PatientsIntraOp_FF.Code);
% PatientsIntraOp_FF.FF_score = zscore(PatientsIntraOp_FF.FF_score);
% 
% % %writetable(PatientsIntraOp_FF,'/home/konstantin/Documents/Data/JR/JR_flipped/FF_Train_Table_V2_low_sigma.csv')
% % 
% % glme = fitglme(PatientsIntraOp_FF,'Response~FF_score+(1|Code)','Distribution','Binomial','Link','logit');
% % 
% let's iterate over test with effect and create a table
I_test = Improvement(test);
counter = 0;
for pt_j = 1:length(gpatsel)

    if test(pt_j) == 1

        counter = counter + 1;

        temp2=strsplit(obj.M.patient.list{gpatsel(pt_j)},'_');
        pt_ID2 = temp2{2}; 
        %trajDepth2 = temp2{end}(4:end);
        side2 = temp2{3};
        side_label = side2;
        
        if contains(obj.M.patient.list{gpatsel(pt_j)}, '_fl.')
            amp = str2num(temp2{end-1}(1:3));
            trajDepth2 = temp2{end-1}(4:end);
        else
            amp = str2num(temp2{end}(1:3));
            trajDepth2 = temp2{end}(4:end);
        end

        row_to_add = {pt_ID2,trajDepth2,side_label,Ihat(pt_j),amp, Improvement(pt_j)};
        if ~exist('PatientsIntraOpTest_FF','var')
            PatientsIntraOpTest_FF = table;
            PatientsIntraOpTest_FF.Code = pt_ID2;
            PatientsIntraOpTest_FF.TrajDepth = trajDepth2;
            PatientsIntraOpTest_FF.Side = side_label;
            PatientsIntraOpTest_FF.FF_score = Ihat(pt_j);
            PatientsIntraOpTest_FF.ClinAny = amp;
            PatientsIntraOpTest_FF.Response = Improvement(pt_j);
        else
            PatientsIntraOpTest_FF = [PatientsIntraOpTest_FF;row_to_add];
        end
    end
end

PatientsIntraOpTest_FF.Code = categorical(PatientsIntraOpTest_FF.Code);
PatientsIntraOpTest_FF.FF_score = zscore(PatientsIntraOpTest_FF.FF_score);
% %writetable(PatientsIntraOpTest_FF,'/home/konstantin/Documents/Data/JR/JR_flipped/FF_Test_Table_V2_low_sigma.csv')
% 
glme2 = fitglme(PatientsIntraOpTest_FF,'Response~FF_score+(1|Code)','Distribution','Binomial','Link','logit','FitMethod','Laplace');

% figure
% plotResiduals(glme2,'fitted','ResidualType','Pearson')
% % 
% predicted_response = predict(glme,PatientsIntraOpTest_FF);
% predicted_response2 = predict(glme2,PatientsIntraOpTest_FF);
% predicted_response3 = predict(mdl_test,PatientsIntraOpTest_FF);
% 
% bin_resp = predicted_response >= scores_thresh;
% bin_resp2 = predicted_response2 >= 0.5;
% bin_resp2 = predicted_response2 >= 0.5;
% false_positives = 0;
% false_negatives = 0;
% false_positives2 = 0;
% false_negatives2 = 0;
% false_positives3 = 0;
% false_negatives3 = 0;
% 
% for i = 1:size(PatientsIntraOpTest_FF,1)
% 
%     if PatientsIntraOpTest_FF.Response(i) == 0 && bin_resp(i) == 1
%         false_positives = false_positives + 1;
%     elseif PatientsIntraOpTest_FF.Response(i) == 1 && bin_resp(i) == 0
%         false_negatives = false_negatives + 1;
%     end
% 
% %     if PatientsIntraOpTest_FF.Response(i) == 0 && bin_resp2(i) == 1
% %         false_positives2 = false_positives2 + 1;
% %     elseif PatientsIntraOpTest_FF.Response(i) == 1 && bin_resp2(i) == 0
% %         false_negatives2 = false_negatives2 + 1;
% %     end
% end
% 
% % first, we fit a logit function for our binary prediction
% mdl_test = fitglm(PatientsIntraOpTest_FF,'Response~FF_score','Distribution','binomial','Link','logit');
% 
% % second, we run ROC curve analysis
% scores = mdl_test.Fitted.Probability;
% scores = glme.Fitted;
% [X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);
% if disp_regression_plots
%     figure, plot(X,Y, 'k', 'linew', 1.5)
%     set(gcf,'color','w');
%     hold on
%     plot(OPTROCPT(1),OPTROCPT(2),'ro', 'MarkerSize',10)
%     xlabel('False positive rate') 
%     ylabel('True positive rate')
%     txt = ['AUC: ' num2str(AUC)];
%     text(0.7,0.1,txt)
%     title('ROC for Classification by Logistic Regression')
% end
% % optimal threshold on the classifier
% scores_thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
% 
% % mdl = fitglm(zscore(Ihat_train),Improvement(training),'Distribution','binomial','Link','logit');

%mdl = fitglm(Ihat_train,Improvement(training),'Distribution','binomial','Link','logit');
mdl = fitglm(Ihat(test),Improvement(test),'Distribution','binomial','Link','logit');

% second, we run ROC curve analysis
%scores = mdl.Fitted.Probability;
%[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);
scores = glme2.Fitted;
[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(test),scores,1);
if disp_regression_plots
    figure, plot(X,Y, 'k', 'linew', 1.5)
    set(gcf,'color','w');
    hold on
    plot(OPTROCPT(1),OPTROCPT(2),'ro', 'MarkerSize',10)
    xlabel('False positive rate') 
    ylabel('True positive rate')
    txt = ['AUC: ' num2str(AUC)];
    text(0.7,0.1,txt)
    title('ROC for Classification by Logistic Regression')
end

% optimal threshold on the classifier
scores_thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));

% plot logit fit for training
vec_val = min(Ihat_train):1:max(Ihat_train);
if disp_regression_plots
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
%scores_test = predict(mdl,Ihat(test));
scores_test = predict(glme2,PatientsIntraOpTest_FF);


Ihat_prediction = scores_test > scores_thresh;

% get the confussion matrix (this can be done on the test set now)

if disp_regression_plots
    figure
    cm = confusionchart(logical(Improvement(test)), Ihat_prediction);
    set(gcf,'color','w');
    
    tp = sum((Ihat_prediction == 1) & (Improvement(test) == 1));
    fp = sum((Ihat_prediction == 1) & (Improvement(test) == 0));
    tn = sum((Ihat_prediction == 0) & (Improvement(test) == 0));
    fn = sum((Ihat_prediction == 0) & (Improvement(test) == 1));
    
    sensitivity = tp/(tp + fn);  % TPR
    specificity = tn/(tn + fp);  % TNR
    
    cm.Title = ['Sensitivity: ', sprintf('%.2f',sensitivity), '; ', 'Specificity: ', sprintf('%.2f',specificity)];
end


% correlation plots for each patient

% vals as contrast of two protocols
%amp_patients = cell(1, 100);  % store amps for each patient
%counter = 1;
I_test = Improvement(test);

pt_i = 1;
difs = -1*ones(100,1);
difs_mean = -1*ones(100,1);
corr_vals_r = zeros(100,1);
corr_vals_p = zeros(100,1);
pt_N = 1;
pt_K = 1;

depths_max = -1*ones(39,1);

lost_prediction = 0;
total_predictions = 0;
amps_not_predicted = [];
all_predicted_amps = [];

while pt_i <= length(gpatsel)-1
    temp=strsplit(obj.M.patient.list{gpatsel(pt_i)},'_');
    %disp(temp{end}(1:3))

    % iterate only over that patients VTAs
    pt_ID = temp{2};
    %trajDepth = temp{end}(4:end);
    %side = temp{3};

    amps = [];
    amps_predicted = [];

    depths = [];  % - will be neglected

    scores = [];
    scores_predicted = [];
    scores_av = [];

    % check both sides
    sides = {'right', 'left'};
    for side_idx = 1:length(sides) 

        true_resp = 0;

        amp_checked = ["None"];
        amp_pred_checked = ["None"];
        amps_side = [];
        amps_predicted_side = [];

        side = sides{side_idx};
        % check amplitude for side-effect
        for pt_j = pt_i:length(gpatsel)
            temp2=strsplit(obj.M.patient.list{gpatsel(pt_j)},'_');
            pt_ID2 = temp2{2}; 
            %trajDepth2 = temp2{end}(4:end);
            side2 = temp2{3};
    
            if contains(obj.M.patient.list{gpatsel(pt_j)}, '_fl.')
                amp = str2num(temp2{end-1}(1:3));
                trajDepth2 = temp2{end-1}(4:end);
            else
                amp = str2num(temp2{end}(1:3));
                trajDepth2 = temp2{end}(4:end);
            end
    
            % exit if different patient
            if ~strcmp(pt_ID2, pt_ID)
                break
            elseif ~strcmp(side,side2)
                continue
            else

                if I_test(pt_j) == 1 && ~any(contains(amp_checked, trajDepth2)) 

                    if strcmp(pt_ID, 'MNI/SBNC46CB')
                        disp("here")
                    end

%                     if amp > 6.0
%                         amp_checked =  [amp_checked,"No response"];
%                         depths = [depths, -1];
%                         amps_side = [amps_side, NaN];
%                     else

    % 
    %                     if amp == 0.5
    %                         disp("____________-")
    %                         disp(pt_ID2)
    %                         disp(trajDepth2)
    %                         disp(side2)
    %                     end
    
                        amps_side = [amps_side, amp];
                        scores = [scores, scores_test(pt_j)];
                        %disp(trajDepth2)
                        if contains(trajDepth2, '.nii')
                            depths = [depths, str2num(trajDepth2(end-6:end-4))];
                        else
                            depths = [depths, str2num(trajDepth2(end-2:end))];
                        end
    
                        %disp(pt_j)
                        %disp(depths(end))
                        %disp("-----------")
    
                        amp_checked =  [amp_checked,trajDepth2];
                        %found_amp = 1;
                        true_resp = true_resp + 1;
                    %end

                elseif amp == 7.0 && I_test(pt_j) == 0 
                    amps_side = [amps_side, NaN];
                    depths = [depths, -1];
                    amp_checked =  [amp_checked,"No response"];
                end
    
                % also skip checking the same elecrode
                if Ihat_prediction(pt_j) == 1 && ~any(contains(amp_pred_checked, trajDepth2))
    
                    amps_predicted_side = [amps_predicted_side, amp];  % predicted amps for response
                    scores_predicted = [scores_predicted, scores_test(pt_j)];
    
                    amp_pred_checked =  [amp_pred_checked,trajDepth2];
                    %found_pred_amp = 1;
    
                    %disp(trajDepth_current2)
                elseif amp == 7.0 && Ihat_prediction(pt_j) == 0
                    amps_predicted_side = [amps_predicted_side, NaN];
                    amp_pred_checked =  [amp_pred_checked,"No prediction"];
                end
            end
        end

        % I need to sort here, only matching pairs
        amps_predicted_sorted_side = [];
        amps_sorted_side = [];
        for i = 2:length(amp_checked)
            if ~isnan(amps_side(i-1))
                for j = 2:length(amp_pred_checked)
                    if ~isnan(amps_predicted_side(j-1)) && strcmp(amp_pred_checked(j), amp_checked(i))
                        amps_sorted_side = [amps_sorted_side, amps_side(i-1)];
                        amps_predicted_sorted_side = [amps_predicted_sorted_side, amps_predicted_side(j-1)];
                    end
                end
            end
        end

        %lost_prediction = lost_prediction + (length(amps_side)-length(amps_sorted_side));
        lost_prediction = lost_prediction + (true_resp-length(amps_sorted_side));
        total_predictions = total_predictions + true_resp;

         if strcmp(pt_ID, 'MNI/SBNC46CB')
             disp("check")
         end

        amps = [amps,amps_sorted_side];
        amps_predicted = [amps_predicted,amps_predicted_sorted_side];
    end


%     amps_mean = round(mean(amps));
% 
%     % check scores at specific amplitude
%     for pt_j = pt_i:length(gpatsel)
%         temp2=strsplit(obj.M.patient.list{gpatsel(pt_j)},'_');
%         pt_ID2 = temp2{2}; 
%         if contains(obj.M.patient.list{gpatsel(pt_j)}, '_fl.')
%             amp= str2num(temp2{end-1}(1:3));
%         else
%             amp = str2num(temp2{end}(1:3));
%         end
% 
%         % exit if different patient
%         if ~strcmp(pt_ID2, pt_ID)
%             break
%         else
%             if amps_mean  == amp
% 
%                 %amps = [amps, amp];
%                 scores_av = [scores_av, scores_test(pt_j)];
%             end
%         end
%     end


    % correlate within the patient
    %figure
    %ea_corrplot(amps,scores)
    %ea_corrplot(amps,amps_predicted,nan, {pt_ID, 'Amp', 'Pred Amp (res 0.5)'})
    %ea_corrplot(amps,scores_av)

    if length(amps(~isnan(amps))) > 3 && length(amps_predicted(~isnan(amps_predicted))) > 3
        [corr_vals_r(pt_N), corr_vals_p(pt_N)] = ea_permcorr(amps',amps_predicted','spearman');
        %[difs(pt_N),idx_max] = ea_nanmax(abs(amps'-amps_predicted'));
        %depths_max(pt_N) = depths(idx_max);
        %difs_mean(pt_N) = ea_nanmean(abs(amps'-amps_predicted'));
        %disp(pt_N)

%         if length(amps_predicted(~isnan(amps_predicted))) > 5 && length(amps(~isnan(amps))) > 5
%             figure
%             ea_corrplot(amps,amps_predicted,nan, {pt_ID, 'Amp', 'Pred Amp (res 0.5)'})
%         end

        if corr_vals_r(pt_N) < 0.33
            disp("here")
        end

        pt_N = pt_N + 1;
    else
        disp("skipped corr")
    end

    %disp(pt_ID)
    if all(isempty(amps))
        difs(pt_K) = [];
        difs_mean(pt_K) = [];
    else
        [difs(pt_K),idx_max] = ea_nanmax(abs(amps'-amps_predicted'));
        difs_mean(pt_K) = ea_nanmean(abs(amps'-amps_predicted'));

        if difs_mean(pt_K) > 2
            disp(pt_ID)
            disp(amps)
        end

    end
    pt_K = pt_K + 1;

    pt_i = pt_j;

    all_predicted_amps = [all_predicted_amps;amps_predicted'];

%     % check only the response ones
%     if Improvement(test(pt_i))  == 1
% 
%         if contains(obj.M.patient.list{gpatsel(pt_i)}, '_fl.')
%             amp{counter} = str2num(temp{end-1}(1:3));
%         else
%             amp{counter} = str2num(temp{end}(1:3));
%         end
%     end

end

%disp("depths_max")
%disp(depths_max)

%figure
%histogram(all_predicted_amps,[0,1,2,3,4,5,6,7,8,9])

disp("Number of lost electrode predictions")
disp(lost_prediction)
disp("Out of ")
disp(total_predictions)

%figure()
%histogram(depths_max,9)

corr_vals_r(pt_N:end) = [];
difs(pt_K:end) = [];
difs_mean(pt_K:end) = [];

% make bar plot of max errors and correlations
small_r =  corr_vals_r(corr_vals_r<0.33);
moderate_r =  corr_vals_r(corr_vals_r>0.33 & corr_vals_r<0.66);
high_r =  corr_vals_r(corr_vals_r>0.66);

X = categorical({'R<0.33', '0.33<R<0.66', 'R>0.66'});
X = reordercats(X,{'R<0.33', '0.33<R<0.66', 'R>0.66'});
Y = [length(small_r), length(moderate_r), length(high_r)];
figure()
bar(X,Y)
title('Fiber Model: Spearman Predicted vs Actual Threshold')


difs_1mA = difs(difs<=1);
difs_2mA = difs(1<difs & difs<=2);
difs_2_plus_mA = difs(2<difs);


X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [length(difs_1mA), length(difs_2mA), length(difs_2_plus_mA)];
figure()
bar(X,Y)
title('Fiber Model: Max Patient Error')

difs_mean_1mA = difs_mean(difs_mean<=1);
difs_mean_2mA = difs_mean(1<difs_mean & difs_mean<=2);
difs_mean_2_plus_mA = difs_mean(2<difs_mean);


X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [length(difs_mean_1mA), length(difs_mean_2mA), length(difs_mean_2_plus_mA)];
figure()
bar(X,Y)
title('Fiber Model: Mean Patient Error')


% % plot for Jan
% X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
% Y = [0, 6, 26 ;length(difs_1mA), length(difs_2mA), length(difs_2_plus_mA)];
% figure()
% bar(X,Y)
% title('Fiber Model: Max in Patient Error')
% legend({'Amplitude Model','Fiber Model'},...
%     'Location','northwest')
% 
% X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
% X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
% Y = [9, 27, 6;length(difs_mean_1mA), length(difs_mean_2mA), length(difs_mean_2_plus_mA)];
% figure()
% bar(X,Y)
% title('Fiber Model: Mean in Patient Error')
% legend('Amplitude Model','Fiber Model')

end
