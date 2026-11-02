function [z,beta,se,p,isBinary] = ea_discfibers_mixed_stat(activ_prob, y, amp, probe_id, varargin)
% Fiber-wise mixed-effects Wald statistic for y ~ x_f + amp + (1|probe)
%   activ_prob : nObs x nFibers
%   y          : nObs x 1, binary (0/1, logical) or continuous
%   amp        : nObs x 1 amplitude
%   probe_id   : nObs x 1 cluster labels (numeric, string, or categorical)
% Name-value: 'Family' = 'auto' (default) | 'binary' | 'continuous'

ip = inputParser;
addParameter(ip,'Family','auto');
parse(ip,varargin{:});

y = y(:); amp = amp(:); probe_id = probe_id(:);
switch lower(ip.Results.Family)
    case 'binary',     isBinary = true;
    case 'continuous', isBinary = false;
    otherwise,         isBinary = islogical(y) || all(ismember(y(~isnan(y)),[0 1]));
end

nF = size(activ_prob,2);
z = nan(1,nF); beta = z; se = z; p = z;

tbl = table(double(y), categorical(probe_id), (amp-mean(amp))/std(amp), ...
            'VariableNames',{'y','probe','amp'});
formula = 'y ~ x + amp + (1|probe)';

parfor f = 1:nF                          % plain 'for' without Parallel Toolbox
    t = tbl;
    t.x = activ_prob(:,f);
    if std(t.x) < 1e-8, continue; end    % fiber never varies -> no information
    ws = warning('off','all');
    try
        if isBinary
            m = fitglme(t, formula, 'Distribution','Binomial', ...
                        'Link','logit', 'FitMethod','Laplace');
        else
            m = fitlme(t, formula);
        end
        c = m.Coefficients;
        r = strcmp(c.Name,'x');
        beta(f) = c.Estimate(r);
        se(f)   = c.SE(r);
        p(f)    = c.pValue(r);
        z(f)    = beta(f)/se(f);
    catch
        % non-convergence or separation: stays NaN
    end
    warning(ws);
end