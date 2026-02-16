function ea_corrplot_ICC(x,y,alpha_vec,xlbl,ylbl)

% 1. Calculate Correlations
% Note: Requires Statistics and Machine Learning Toolbox
[rho_p, p_p] = corr(x, y, 'type', 'Pearson');
[rho_s, p_s] = corr(x, y, 'type', 'Spearman');

baseColor = [0 0.447 0.741];

% 2. Create the Plot
figure('Color', 'w');
hold on;

% Loop to apply the alpha vector to each point individually
for i = 1:length(x)
    scatter(x(i), y(i), 200, 'MarkerFaceColor', [0 0.447 0.741], ...
            'MarkerEdgeColor', baseColor, 'MarkerFaceAlpha', alpha_vec(i), ...
            'LineWidth', 1.5);
end

% 3. Add a basic regression line for visual correlation
coeffs = polyfit(x, y, 1);
x_fit = linspace(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 100);
y_fit = polyval(coeffs, x_fit);
plot(x_fit, y_fit, 'r--', 'LineWidth', 1.5);

customMap = [linspace(1, baseColor(1), 256)', ... 
    linspace(1, baseColor(2), 256)', ... 
    linspace(1, baseColor(3), 256)']; 
colormap(customMap); cb = colorbar; 
ylabel(cb, 'ICC weight'); 
caxis([0 1]); % Force scale from 0 to 1

% 4. Format the plot
xlabel(xlbl);
ylabel(ylbl);
% 5. Add Correlation Results as a Text Box
% stats_str = {sprintf('Pearson r: %.3f', rho_p), ...
%              sprintf('Spearman \\rho: %.3f', rho_s)};

stats_str = sprintf('Pearson r=%.3f p=%.5f   Spearman \\rho: %.3f p=%.5f', rho_p,p_p, rho_s,p_s);


title(stats_str);
grid on;
set(gca, 'FontSize', 9);

         
% text(0.05, 0.95, stats_str, 'Units', 'normalized', ...
%     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ... 
%     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'latex', ... 
%     'FontSize', 12);

hold off;