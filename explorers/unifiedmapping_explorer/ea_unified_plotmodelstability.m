function h = ea_unified_plotmodelstability(S, parentHandle)
% ea_unified_plotmodelstability  Plot ea_unified_modelstability.m output.
%
% h = ea_unified_plotmodelstability(S)               -- opens a standalone figure
% h = ea_unified_plotmodelstability(S, parentHandle) -- parents into an existing
%                                                uipanel/uifigure instead
%
% A 2x4 grid of seven panels (last tile left empty): top row compares each
% model to the PREVIOUS sample size (Spearman, Dice, centroid
% displacement); bottom row compares each model directly to the model fit
% on ALL patients instead (Spearman/Dice/centroid vs-final), plus
% nSelected. Mean line across repeats with a +/-SD ribbon, x = sample size,
% y-axis labelled on every panel. Panels with a defined stability
% criterion also mark the +/-tol band and the stability point. Centroid
% panels are replaced with a placeholder if S.centroid is empty (no
% coordinates were supplied).
%
% Returns a struct of handles (figure/tiledlayout/axes/lines/ribbons) so
% the caller can restyle.

if nargin < 2 || isempty(parentHandle)
    parentHandle = figure('Color','w', 'Name','Model Stability', ...
        'NumberTitle','off', 'Position',[100 100 1800 650]);
end

col       = [0.20 0.45 0.70]; % single sequential hue for every mean line
gridCol   = [0.85 0.85 0.85];
axCol     = [0.40 0.40 0.40];
bandCol   = [0.75 0.75 0.75];
lineCol   = [0.55 0.55 0.55];
labelCol  = [0.45 0.45 0.45];

nPanels = 7;
tl = tiledlayout(parentHandle, 2, 4, 'Padding','compact', 'TileSpacing','compact');

h = struct;
h.figureOrParent = parentHandle;
h.tiledlayout = tl;
h.ax     = gobjects(1,nPanels);
h.mean   = gobjects(1,nPanels);
h.ribbon = gobjects(1,nPanels);

% {stabilityPoint field name, mean curve, sd curve, title, ylabel}
panelSpecs = { ...
    'spearman',            S.spearman_mean,             S.spearman_sd,             'Dense-vector Spearman (successive)',      'Spearman \rho'; ...
    'dice',                S.dice_mean.pooled,          S.dice_sd.pooled,          'Selected-set Dice (successive)',          'Dice coefficient'; ...
    'centroid',            S.centroid_mean,             S.centroid_sd,             'Centroid displacement (successive)',      'Displacement (mm)'; ...
    'nSelected',           S.nSelected_mean,            S.nSelected_sd,            '# elements selected',                     '# elements'; ...
    'spearman_vs_final',   S.spearman_vs_final_mean,    S.spearman_vs_final_sd,    'Dense-vector Spearman (vs. final model)', 'Spearman \rho'; ...
    'dice_vs_final',       S.dice_vs_final_mean.pooled, S.dice_vs_final_sd.pooled, 'Selected-set Dice (vs. final model)',     'Dice coefficient'; ...
    'centroid_vs_final',   S.centroid_vs_final_mean,    S.centroid_vs_final_sd,    'Centroid displacement (vs. final model)', 'Displacement (mm)' ...
};

x = S.sizes(:);

for p = 1:nPanels
    ax = nexttile(tl);
    hold(ax, 'on');
    h.ax(p) = ax;

    name = panelSpecs{p,1};
    m    = panelSpecs{p,2};
    s    = panelSpecs{p,3};
    ttl  = panelSpecs{p,4};
    ylab = panelSpecs{p,5};

    if isempty(m)
        text(ax, 0.5, 0.5, 'No coordinates supplied', ...
            'HorizontalAlignment','center', 'Units','normalized', 'Color', labelCol);
        axis(ax, 'off');
        title(ax, ttl, 'Color', [0.25 0.25 0.25], 'FontWeight','normal');
        continue
    end

    valid = ~isnan(m);
    xv = x(valid); mv = m(valid); sv = s(valid);
    sv(isnan(sv)) = 0;

    if numel(xv) >= 2
        h.ribbon(p) = fill(ax, [xv; flipud(xv)], [mv-sv; flipud(mv+sv)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    h.mean(p) = plot(ax, xv, mv, '-o', 'Color', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1.5, 'MarkerSize', 4);

    if isfield(S.stabilityPoint, name) && ~isnan(S.stabilityPoint.(name)) && ~isempty(mv)
        sp = S.stabilityPoint.(name);
        finalVal = mv(end);
        yline(ax, finalVal*(1-S.opts.tol), ':', 'Color', bandCol);
        yline(ax, finalVal*(1+S.opts.tol), ':', 'Color', bandCol);
        xline(ax, sp, '--', 'Color', lineCol);
        yl = ylim(ax); xl = xlim(ax);
        if sp > mean(xl)
            halign = 'right'; labelStr = sprintf('stable @ n=%d ', sp);
        else
            halign = 'left'; labelStr = sprintf(' stable @ n=%d', sp);
        end
        text(ax, sp, yl(2), labelStr, 'Color', labelCol, ...
            'FontSize', 8, 'VerticalAlignment','top', 'HorizontalAlignment', halign);
    end

    ax.XColor = axCol; ax.YColor = axCol;
    ax.GridColor = gridCol; grid(ax, 'on'); box(ax,'off');
    xlabel(ax, 'Sample size (n patients)');
    ylabel(ax, ylab);
    title(ax, ttl, 'Color', [0.25 0.25 0.25], 'FontWeight','normal');
end

end
