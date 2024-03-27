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
        %side2 = temp2{3};
        %side_label = side2;
        
        if contains(obj.M.patient.list{gpatsel(pt_j)}, '_fl.')
            %amp = str2num(temp2{end-1});
            %trajDepth2 = temp2{end}(4:end);
            side2 = 'left';
        else
            %amp = str2num(temp2{end});
            %trajDepth2 = temp2{end}(4:end);
            side2 = 'right';
        end

        trajDepth2 = temp2{end}(4:end);
        side_label = side2;
        amp = str2num(temp2{end-1});

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
%PatientsIntraOpTest_FF.FF_score = zscore(PatientsIntraOpTest_FF.FF_score);
% %writetable(PatientsIntraOpTest_FF,'/home/konstantin/Documents/Data/JR/JR_flipped/FF_Test_Table_V2_low_sigma.csv')
% 
glme = fitglme(PatientsIntraOpTest_FF,'Response~FF_score+(1|Code)','Distribution','Binomial','Link','logit','FitMethod','Laplace');



% % let's predict post-op
% 
% counter = 0;
% for pt_j = 1:length(gpatsel)
% 
%     if test(pt_j) == 1
% 
%         counter = counter + 1;
% 
%         temp2=strsplit(obj.M.patient.list{gpatsel(pt_j)},'_');
%         pt_ID2 = temp2{2}; 
%         %trajDepth2 = temp2{end}(4:end);
%         if contains(obj.M.patient.list{gpatsel(pt_j)},'left')
%             side2 = 'left';
%         else
%             side2 = 'right';
%         end
%         %side2 = temp2{3};
%         side_label = side2;
%         
%         amp = temp2{end-1};
%         contact_N = temp2{end}(1);
% 
%         row_to_add = {pt_ID2,contact_N,side_label,Ihat(pt_j),amp};
%         if ~exist('PatientsPostOpTest_FF','var')
%             PatientsPostOpTest_FF = table;
%             PatientsPostOpTest_FF.Code = pt_ID2;
%             PatientsPostOpTest_FF.Contact = contact_N;
%             PatientsPostOpTest_FF.Side = side_label;
%             PatientsPostOpTest_FF.FF_score = Ihat(pt_j);
%             PatientsPostOpTest_FF.Amp = amp;
%             %PatientsPostOpTest_FF.Response = Improvement(pt_j);
%         else
%             PatientsPostOpTest_FF = [PatientsPostOpTest_FF;row_to_add];
%         end
%     end
% end
% 
% PatientsPostOpTest_FF.Code = categorical(PatientsPostOpTest_FF.Code);
% %PatientsPostOpTest_FF.FF_score = zscore(PatientsPostOpTest_FF.FF_score);
% scores_test = predict(glme2,PatientsPostOpTest_FF);
% Ihat_prediction = scores_test > 0.5913;
% %Ihat_prediction = scores_test > 0.7023;
% PatientsPostOpTest_FF.Response = Ihat_prediction;

%scores_test = predict(mdl,Ihat(test));


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
%[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);
%scores = glme.Fitted;
scores = mdl.Fitted.Probability;
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
scores_test = predict(mdl,Ihat(test));
%scores_test = predict(glme,PatientsIntraOpTest_FF);


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



% let's check prediction accuracy for each patient
predicted_response = Ihat_prediction;

% let's find the actual amplitudes
% I will need to re-add again from the original file checking for STN
load('/home/konstantin/Documents/Data/JR/PatientsEMGOR_corr.mat')


pt_counter = 0;
pt_label = 'JohnDoe';
difs = cell(42,1);
difs_1mA = 0;
difs_2mA = 0;
difs_2_plus_mA = 0;
difs_mean_1mA = 0;
difs_mean_2mA = 0;
difs_mean_2_plus_mA = 0;
already_checked = [];
for measure_i = 1:length(predicted_response)

    % only check if over the threshold and response positive
    if predicted_response(measure_i) == 1

        parts = strsplit(obj.M.patient.list{gpatsel(measure_i)},'_');
        % stupid way to recover label and traj depth
        trajDepth = parts{end}(4:end-4);
        pt_checked = parts{2}(5:end);
        amp = str2num(parts{end-1});
    
        if contains(obj.M.patient.list{gpatsel(measure_i)}, 'right')
            side_label = 'right';
            hemi_label = 're';
        else
            side_label = 'left';
            hemi_label = 'li';
        end
    
        traj = char(extract(trajDepth,lettersPattern));
        depth = str2num(trajDepth(length(traj)+1:end));
    
        % let's find clinical thresholds for this patient, trajectory and side
        pt_entry = find(PatientsEMGOR4.Code == pt_checked & PatientsEMGOR4.Hemisphere == hemi_label & PatientsEMGOR4.Traj == traj & PatientsEMGOR4.Depth == depth);
        
        if ismember(pt_entry,already_checked)
            %already checked this electrode
            continue
        else
            already_checked = [already_checked,pt_entry];
        end

        if length(pt_entry) > 1
            disp(vta_i)
        end
    
        if ~isinf(PatientsEMGOR4.ClinAny(pt_entry)) && ~isnan(PatientsEMGOR4.ClinAny(pt_entry))
            amp_actual = PatientsEMGOR4.ClinAny(pt_entry);
        end


        if ~strcmp(pt_checked,pt_label)
            % new patient
            pt_label = pt_checked;
            pt_counter = pt_counter + 1;
            difs{pt_counter} = abs(amp_actual - amp);
        else
            difs{pt_counter} = [difs{pt_counter},abs(amp_actual - amp)];
        end
    end
end

% you can check how many indices are missing in already_checked

for pt_i = 1:42
    if ~isempty(find(2<difs{pt_i}))
        difs_2_plus_mA = difs_2_plus_mA + 1;
    elseif ~isempty(find(1<difs{pt_i} & difs{pt_i}<=2))
        difs_2mA = difs_2mA + 1;
    else
        difs_1mA = difs_1mA + 1;
    end
    if 2<mean(difs{pt_i})
        difs_mean_2_plus_mA = difs_mean_2_plus_mA + 1;
    elseif 1<mean(difs{pt_i}) && mean(difs{pt_i})<=2
        difs_mean_2mA = difs_mean_2mA + 1;
    else
        difs_mean_1mA = difs_mean_1mA + 1;
    end
    % difs_1mA = difs_1mA + length(find(difs{pt_i}<=1));
    % difs_2mA = difs_2mA + length(find(1<difs{pt_i} & difs{pt_i}<=2));
    % difs_2_plus_mA = difs_2_plus_mA + length(find(2<difs{pt_i}));
end
X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [difs_1mA, difs_2mA, difs_2_plus_mA];
figure()
bar(X,Y)
title('Fiber Model: Max Patient Error')
X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [difs_mean_1mA, difs_mean_2mA, difs_mean_2_plus_mA];
figure()
bar(X,Y)
title('Fiber Model: Mean Patient Error')

end
