function [z_score, p_value] = delong_test_independent(labels1, scores1, labels2, scores2)
    % Calculates DeLong test for two independent AUCs
    
    % 1. Calculate AUCs and Components for Model 1
    [auc1, V10, V01] = calculate_delong_components(labels1, scores1);
    var1 = (var(V10)/length(V10)) + (var(V01)/length(V01));
    
    % 2. Calculate AUCs and Components for Model 2
    [auc2, V10_2, V01_2] = calculate_delong_components(labels2, scores2);
    var2 = (var(V10_2)/length(V10_2)) + (var(V01_2)/length(V01_2));
    
    % 3. Compare (Standard Error for independent samples)
    se_diff = sqrt(var1 + var2);
    z_score = (auc1 - auc2) / se_diff;
    p_value = 2 * (1 - normcdf(abs(z_score))); % Two-tailed test
    
    fprintf('AUC 1: %.4f\nAUC 2: %.4f\n', auc1, auc2);
    fprintf('Z-score: %.4f\nP-value: %.4f\n', z_score, p_value);
end

function [auc, V10, V01] = calculate_delong_components(labels, scores)
    % Get logical indices
    pos = scores(labels == 1);
    neg = scores(labels == 0);
    m = length(pos);
    n = length(neg);
    
    % Structural components
    V10 = zeros(m, 1);
    V01 = zeros(n, 1);
    
    for i = 1:m
        V10(i) = sum(neg < pos(i)) / n + 0.5 * sum(neg == pos(i)) / n;
    end
    for j = 1:n
        V01(j) = sum(pos > neg(j)) / m + 0.5 * sum(pos == neg(j)) / m;
    end
    
    auc = mean(V10);
end