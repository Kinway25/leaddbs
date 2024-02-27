function ea_get_binary_prediction_plots(obj, N_resp, gpatsel, Ihat_prediction, I_test)

pt_i = 1;
difs = -1*ones(100,1);
difs_mean = -1*ones(100,1);
corr_vals_r = zeros(100,1);
corr_vals_p = zeros(100,1);
pt_N = 1;
pt_K = 1;

lost_prediction = 0;

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

% 
%                     if amp == 0.5
%                         disp("____________-")
%                         disp(pt_ID2)
%                         disp(trajDepth2)
%                         disp(side2)
%                     end

                    amps_side = [amps_side, amp];
                    %scores = [scores, scores_test(pt_j)];
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
                elseif amp == 7.0 && I_test(pt_j) == 0
                    amps_side = [amps_side, NaN];
                    depths = [depths, -1];
                    amp_checked =  [amp_checked,"No response"];
                end
    
                % also skip checking the same elecrode
                if Ihat_prediction(pt_j) == 1 && ~any(contains(amp_pred_checked, trajDepth2))
    
                    amps_predicted_side = [amps_predicted_side, amp];  % predicted amps for response
                    %scores_predicted = [scores_predicted, scores_test(pt_j)];
    
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

        lost_prediction = lost_prediction + (length(amps_side)-length(amps_sorted_side));

%         if strcmp(pt_ID, 'MNI/sMKbfUxmYtXmA2q6KxRveR')
%             disp("check")
%         end

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

        if corr_vals_r(pt_N) < 0.33
            disp("here")
        end

        pt_N = pt_N + 1;
    else
        disp(pt_ID)
    end

    %disp(pt_ID)
    [difs(pt_K),idx_max] = ea_nanmax(abs(amps'-amps_predicted'));
    difs_mean(pt_K) = ea_nanmean(abs(amps'-amps_predicted'));
    pt_K = pt_K + 1;

    pt_i = pt_j;

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

disp("Number of lost electrode predictions")
disp(lost_prediction)

%figure()
%histogram(depths_max,9)

corr_vals_r(pt_N:end) = [];
difs(pt_K:end) = [];
difs_mean(pt_K:end) = [];

% make bar plot of max errors and correlations
small_r =  corr_vals_r(corr_vals_r<0.33);
moderate_r =  corr_vals_r(corr_vals_r>0.33 & corr_vals_r<0.66);
high_r =  corr_vals_r(corr_vals_r>0.66);

%disp("depths_max")

%figure()
%histogram(depths_max,9)

corr_vals_r(pt_N:end) = [];
difs(pt_N:end) = [];
difs_mean(pt_N:end) = [];

% make bar plot of max errors and correlations
small_r =  corr_vals_r(corr_vals_r<0.33);
moderate_r =  corr_vals_r(corr_vals_r>0.33 & corr_vals_r<0.66);
high_r =  corr_vals_r(corr_vals_r>0.66);

X = categorical({'Small R', 'Moderate R', 'High R'});
X = reordercats(X,{'Small R', 'Moderate R', 'High R'});
Y = [length(small_r), length(moderate_r), length(high_r)];
figure()
bar(X,Y)
title('Ampl. Model: Spearman Predicted vs Actual Threshold')


difs_1mA = difs(difs<=1);
difs_2mA = difs(1<difs & difs<=2);
difs_2_plus_mA = difs(2<difs);


X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [length(difs_1mA), length(difs_2mA), length(difs_2_plus_mA)];
figure()
bar(X,Y)
title('Ampl. Model: Max Patient Error')


difs_mean_1mA = difs_mean(difs_mean<=1);
difs_mean_2mA = difs_mean(1<difs_mean & difs_mean<=2);
difs_mean_2_plus_mA = difs_mean(2<difs_mean);


X = categorical({'<=1 mA', '<=2 mA', '>2 mA'});
X = reordercats(X,{'<=1 mA', '<=2 mA', '>2 mA'});
Y = [length(difs_mean_1mA), length(difs_mean_2mA), length(difs_mean_2_plus_mA)];
figure()
bar(X,Y)
title('Ampl. Model: Mean Patient Error')