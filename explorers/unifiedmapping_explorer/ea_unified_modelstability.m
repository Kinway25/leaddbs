function S = ea_unified_modelstability(X, y, opts)
% ea_unified_modelstability  Sample-size stability analysis for a mass-univariate
% Unified Mapping Explorer model (fiber filtering or sweet-spot mapping).
%
% Implements the additive-subsample stability curves described in Nowacki
% et al. (Frontiers in Computational Neuroscience, PMC12868180): repeatedly
% grow a random patient subsample and ask how much the fitted per-element
% model changes each time two more patients are added. Also runs a cheap
% leave-one-out (LOO) tier for free, using the n LOO models to (a) estimate
% agreement at n-1 with no extra fitting, and (b) flag individual
% high-leverage patients.
%
% This function is agnostic to whether X holds per-fiber overlap values or
% per-voxel E-field/statistic values -- it only ever sees a dense
% [nElements x nPatients] matrix and an outcome vector.
%
% INPUTS
%   X    - [nElements x nPatients] double. Per-element, per-patient values
%          (fiber overlap for fiber filtering; voxel values for sweet-spot
%          mapping).
%   y    - [nPatients x 1] outcome vector.
%   opts - struct (all fields optional):
%     .statFcn      function handle @(Xsub,ysub) -> [nElements x 1] dense
%                    per-element statistic. Default: vectorized per-row
%                    correlation (Pearson or Spearman, see .corrType).
%                    NOTE: this is a 2-arg/1-output contract, which does
%                    NOT match the 3-arg/2-output signature used by the
%                    ea_explorer_stats_*.m test functions elsewhere in
%                    Lead-DBS ([valsout,psout]=fcn(valsin,outcomein,H0)).
%                    To reuse one of those as .statFcn, wrap it, e.g.:
%                       H0 = obj.statsettings.H0;
%                       opts.statFcn = @(X,y) ea_unified_modelstability_firstout( ...
%                           @() ea_explorer_stats_ranksumtest(X,y,H0));
%                    (ea_unified_modelstability_firstout is a tiny local helper
%                    below, exposed only for copy-paste convenience -- it
%                    is not called anywhere internally.)
%     .coords       [nElements x 3] element coordinates in mm, optional.
%                    Enables the centroid-displacement metric. For fiber
%                    filtering, reduce each streamline in fibcell to one
%                    point first, e.g. mean(fibcell{g,s}{i},1). For
%                    sweet-spot voxels, convert linear mask indices via
%                    ind2sub + ea_vox2mm. Skipped (S.centroid = []) if
%                    omitted.
%     .selection    'topk' (default) or 'percentile'. Selection rule used
%                    UNIFORMLY on every model this function compares, for
%                    the Dice/centroid metrics only -- independent of
%                    whatever threshold the caller's main analysis used.
%     .k            elements per tail for 'topk'. Default 1000.
%     .percentile   percent of covered elements per tail for 'percentile'.
%                    Default 5.
%     .tails        'positive' | 'negative' | 'both'. Default 'both'.
%     .minCoverage  an element is included only if at least this fraction
%                    of patients IN THE CURRENT SUBSAMPLE have a non-zero,
%                    non-NaN value for it. Recomputed within every
%                    subsample (not once globally). Default 0.2.
%     .sizes        subsample sizes to evaluate. Default: additive
%                    schedule 4:2:nPatients (padded with nPatients if the
%                    step doesn't land on it exactly).
%     .nRep         random repeats per size. Default 10.
%     .corrType     'Spearman' (default) or 'Pearson'. Only used by the
%                    default .statFcn.
%     .seed         RNG seed. Default 42.
%     .tol          stability-band tolerance (fraction of final value).
%                    Default 0.05.
%
% OUTPUT S
%   .sizes            [1 x nSizes] sample sizes evaluated.
%   .spearman         [nSizes x nRep] dense-statistic-vector Spearman
%                      agreement between successive models. Row 1 is NaN
%                      (no previous model to compare the first size to).
%                      This metric never uses .selection, so it is
%                      invariant to how the caller's main analysis (or this
%                      function's own Dice/centroid) thresholds anything --
%                      treat it as the primary metric.
%   .spearman_mean/_sd  mean/SD across repeats, omitting NaNs.
%   .dice             struct with .pos/.neg/.pooled, each [nSizes x nRep].
%                      Dice coefficient of the selected (per .selection)
%                      element sets between successive models. .pooled is
%                      the headline number (union of tails).
%   .dice_mean/_sd      same struct shape, mean/SD across repeats.
%   .centroid         [nSizes x nRep] centroid displacement in mm between
%                      successive models' pooled selected sets. Empty ([])
%                      if .coords was not supplied.
%   .centroid_mean/_sd  as above, empty if no coords.
%   .nSelected        [nSizes x nRep] number of pooled selected elements
%                      per model (diagnostic -- shows how set size itself
%                      drifts with n, independent of agreement).
%   .nSelected_mean/_sd
%   .spearman_vs_final  [nSizes x nRep] dense-statistic-vector Spearman
%                      agreement between the model at THIS sample size and
%                      the model fit on ALL patients (.finalModel below) --
%                      as opposed to .spearman, which only ever compares
%                      consecutive sizes. No leading NaN: every valid size
%                      (including the first) has a defined value here,
%                      since the final model is a fixed reference, not a
%                      "previous" model in the growth path. Answers "how
%                      close is this sample size to the end result",
%                      rather than "how much did the last increment
%                      change things".
%   .spearman_vs_final_mean/_sd
%   .dice_vs_final    struct with .pos/.neg/.pooled, each [nSizes x nRep].
%                      Same idea as .dice, but against .finalModel's
%                      selected sets instead of the previous size's.
%   .dice_vs_final_mean/_sd
%   .centroid_vs_final  [nSizes x nRep] centroid displacement in mm between
%                      this size's pooled selected set and .finalModel's.
%                      Empty ([]) if .coords was not supplied.
%   .centroid_vs_final_mean/_sd
%   .stabilityPoint   struct .spearman/.dice/.centroid/.spearman_vs_final/
%                      .dice_vs_final/.centroid_vs_final: earliest sample
%                      size after which that mean curve never leaves a
%                      +/-.tol band around its own final value. NaN if
%                      never reached (or if coords weren't supplied, for
%                      .centroid/.centroid_vs_final).
%   .finalModel       the model fit once on ALL nPatients (not a
%                      subsample): .statVec/.maskIdx (dense statistic and
%                      which elements were covered) and .selPos/.selNeg/
%                      .selPool (the selected element indices under
%                      .selection/.k or .percentile/.tails). Useful when
%                      the caller wants the actual final selected set, not
%                      just its stability across sample sizes.
%   .loo              free LOO tier (see ea_unified_modelstability_loo below):
%                      .pairwise.spearman_mean/_sd, .dice_mean/_sd,
%                      .centroid_mean/_sd (agreement among all n LOO
%                      models, i.e. stability at n-1 with no extra
%                      fitting), and .influence.patientIdx/.score (per-
%                      patient departure of the LOO-without-them model
%                      from the full-sample model, 1-Spearman, sorted
%                      descending -- higher = more influential/high-
%                      leverage patient).
%   .nDroppedSubsamples  count of (size,repeat) combinations skipped
%                      because the subsample was degenerate (zero
%                      coverage, or zero variance in y).
%   .opts             resolved options struct (provenance).
%
% EXAMPLE -- fiber filtering
%   connid   = ea_conn2connid(explorer.calcsettings.fibfilt_connectome);
%   methodid = ea_unifiedmapping_method2methodid(explorer);
%   X = full(explorer.results.fiberfiltering.(connid).(methodid).fibsval{1});
%   X = X(:, explorer.patientselection);
%   y = explorer.responsevar(explorer.patientselection, 1);
%   fibcell = explorer.results.fiberfiltering.(connid).(methodid).fibcell{1};
%   opts = struct('coords', cell2mat(cellfun(@(f) mean(f,1), fibcell, 'Uni', 0)));
%   S = ea_unified_modelstability(X, y, opts);
%
% EXAMPLE -- sweet-spot mapping (voxel-wise, no coordinates)
%   X = explorer.results.sweetspotmapping.efield{1}';
%   X = X(:, explorer.patientselection);
%   y = explorer.responsevar(explorer.patientselection, 1);
%   S = ea_unified_modelstability(X, y);   % opts omitted -> no centroid metric

if nargin < 3 || isempty(opts)
    opts = struct;
end

nPatients = size(X,2);

if size(y,1) ~= nPatients
    ea_error('ea_unified_modelstability: y must have one entry per column of X (nPatients).');
end
y = y(:);

if nPatients < 6
    ea_error('ea_unified_modelstability: need at least 6 patients (nPatients=%d given) -- successive-model comparison is not meaningful below this.', nPatients);
end

opts = ea_unified_modelstability_resolveopts(opts, nPatients);

if numel(opts.sizes) < 2
    ea_error('ea_unified_modelstability: opts.sizes must contain at least 2 sample sizes.');
end

rng(opts.seed);

nSizes = numel(opts.sizes);
nRep   = opts.nRep;
haveCoords = ~isempty(opts.coords);

% Fit once on ALL patients up front so every subsample below can be
% compared directly against it (.spearman_vs_final/.dice_vs_final), not
% just against the previous size.
[finalStat, finalMaskIdx, finalOk] = ea_unified_modelstability_fitone(X, y, opts);
if finalOk
    [finalSelPos, finalSelNeg, finalSelPool] = ea_unified_modelstability_select(finalStat, finalMaskIdx, opts);
else
    finalSelPos = []; finalSelNeg = []; finalSelPool = [];
end

spearman   = nan(nSizes, nRep);
dicePos    = nan(nSizes, nRep);
diceNeg    = nan(nSizes, nRep);
dicePooled = nan(nSizes, nRep);
centroid   = nan(nSizes, nRep);
nSelected  = nan(nSizes, nRep);

spearmanVsFinal   = nan(nSizes, nRep);
dicePosVsFinal    = nan(nSizes, nRep);
diceNegVsFinal    = nan(nSizes, nRep);
dicePooledVsFinal = nan(nSizes, nRep);
centroidVsFinal   = nan(nSizes, nRep);

nDropped = 0;

for r = 1:nRep
    order = randperm(nPatients);

    prevValid   = false;
    prevStat    = [];
    prevMaskIdx = [];
    prevSelPos  = [];
    prevSelNeg  = [];
    prevSelPool = [];

    for i = 1:nSizes
        n_i  = opts.sizes(i);
        idx  = order(1:n_i);
        Xsub = X(:, idx);
        ysub = y(idx);

        [statVec, maskIdx, ok] = ea_unified_modelstability_fitone(Xsub, ysub, opts);

        if ~ok
            nDropped = nDropped + 1;
            continue % this (i,r) stays NaN everywhere; prev* still points to the last valid model
        end

        [selPos, selNeg, selPool] = ea_unified_modelstability_select(statVec, maskIdx, opts);
        nSelected(i,r) = numel(selPool);

        if finalOk
            spearmanVsFinal(i,r) = ea_unified_modelstability_spearman(statVec, maskIdx, finalStat, finalMaskIdx);
            [dicePosVsFinal(i,r), diceNegVsFinal(i,r), dicePooledVsFinal(i,r)] = ea_unified_modelstability_dice( ...
                selPos, selNeg, selPool, finalSelPos, finalSelNeg, finalSelPool);
            if haveCoords
                centroidVsFinal(i,r) = ea_unified_modelstability_centroiddist(opts.coords, selPool, finalSelPool);
            end
        end

        if prevValid
            spearman(i,r) = ea_unified_modelstability_spearman(prevStat, prevMaskIdx, statVec, maskIdx);
            [dicePos(i,r), diceNeg(i,r), dicePooled(i,r)] = ea_unified_modelstability_dice( ...
                prevSelPos, prevSelNeg, prevSelPool, selPos, selNeg, selPool);
            if haveCoords
                centroid(i,r) = ea_unified_modelstability_centroiddist(opts.coords, prevSelPool, selPool);
            end
        end

        prevValid = true;
        prevStat = statVec; prevMaskIdx = maskIdx;
        prevSelPos = selPos; prevSelNeg = selNeg; prevSelPool = selPool;
    end
end

S = struct;
S.sizes = opts.sizes(:)';

S.spearman      = spearman;
S.spearman_mean = mean(spearman, 2, 'omitnan');
S.spearman_sd   = std(spearman, 0, 2, 'omitnan');

S.dice.pos    = dicePos;
S.dice.neg    = diceNeg;
S.dice.pooled = dicePooled;
S.dice_mean.pos    = mean(dicePos, 2, 'omitnan');
S.dice_mean.neg    = mean(diceNeg, 2, 'omitnan');
S.dice_mean.pooled = mean(dicePooled, 2, 'omitnan');
S.dice_sd.pos    = std(dicePos, 0, 2, 'omitnan');
S.dice_sd.neg    = std(diceNeg, 0, 2, 'omitnan');
S.dice_sd.pooled = std(dicePooled, 0, 2, 'omitnan');

if haveCoords
    S.centroid      = centroid;
    S.centroid_mean = mean(centroid, 2, 'omitnan');
    S.centroid_sd   = std(centroid, 0, 2, 'omitnan');
else
    S.centroid      = [];
    S.centroid_mean = [];
    S.centroid_sd   = [];
end

S.nSelected      = nSelected;
S.nSelected_mean = mean(nSelected, 2, 'omitnan');
S.nSelected_sd   = std(nSelected, 0, 2, 'omitnan');

S.spearman_vs_final      = spearmanVsFinal;
S.spearman_vs_final_mean = mean(spearmanVsFinal, 2, 'omitnan');
S.spearman_vs_final_sd   = std(spearmanVsFinal, 0, 2, 'omitnan');

S.dice_vs_final.pos    = dicePosVsFinal;
S.dice_vs_final.neg    = diceNegVsFinal;
S.dice_vs_final.pooled = dicePooledVsFinal;
S.dice_vs_final_mean.pos    = mean(dicePosVsFinal, 2, 'omitnan');
S.dice_vs_final_mean.neg    = mean(diceNegVsFinal, 2, 'omitnan');
S.dice_vs_final_mean.pooled = mean(dicePooledVsFinal, 2, 'omitnan');
S.dice_vs_final_sd.pos    = std(dicePosVsFinal, 0, 2, 'omitnan');
S.dice_vs_final_sd.neg    = std(diceNegVsFinal, 0, 2, 'omitnan');
S.dice_vs_final_sd.pooled = std(dicePooledVsFinal, 0, 2, 'omitnan');

if haveCoords
    S.centroid_vs_final      = centroidVsFinal;
    S.centroid_vs_final_mean = mean(centroidVsFinal, 2, 'omitnan');
    S.centroid_vs_final_sd   = std(centroidVsFinal, 0, 2, 'omitnan');
else
    S.centroid_vs_final      = [];
    S.centroid_vs_final_mean = [];
    S.centroid_vs_final_sd   = [];
end

S.stabilityPoint = struct;
S.stabilityPoint.spearman = ea_unified_modelstability_stabpoint(S.sizes, S.spearman_mean, opts.tol);
S.stabilityPoint.dice     = ea_unified_modelstability_stabpoint(S.sizes, S.dice_mean.pooled, opts.tol);
S.stabilityPoint.spearman_vs_final = ea_unified_modelstability_stabpoint(S.sizes, S.spearman_vs_final_mean, opts.tol);
S.stabilityPoint.dice_vs_final     = ea_unified_modelstability_stabpoint(S.sizes, S.dice_vs_final_mean.pooled, opts.tol);
if haveCoords
    S.stabilityPoint.centroid = ea_unified_modelstability_stabpoint(S.sizes, S.centroid_mean, opts.tol);
    S.stabilityPoint.centroid_vs_final = ea_unified_modelstability_stabpoint(S.sizes, S.centroid_vs_final_mean, opts.tol);
else
    S.stabilityPoint.centroid = nan;
    S.stabilityPoint.centroid_vs_final = nan;
end

S.finalModel = struct('statVec', finalStat, 'maskIdx', finalMaskIdx, 'ok', finalOk, ...
    'selPos', finalSelPos, 'selNeg', finalSelNeg, 'selPool', finalSelPool);

S.loo = ea_unified_modelstability_loo(X, y, opts);

S.nDroppedSubsamples = nDropped;
S.opts = opts;

fprintf('ea_unified_modelstability: %d/%d (size,repeat) subsamples were degenerate and dropped.\n', ...
    nDropped, nSizes*nRep);

end

% ======================================================================
function opts = ea_unified_modelstability_resolveopts(opts, nPatients)
if ~isfield(opts,'statFcn'),     opts.statFcn = [];        end
if ~isfield(opts,'coords') || isempty(opts.coords), opts.coords = []; end
if ~isfield(opts,'selection') || isempty(opts.selection), opts.selection = 'topk'; end
if ~isfield(opts,'k') || isempty(opts.k),                 opts.k = 1000;          end
if ~isfield(opts,'percentile') || isempty(opts.percentile), opts.percentile = 5;  end
if ~isfield(opts,'tails') || isempty(opts.tails),         opts.tails = 'both';    end
if ~isfield(opts,'minCoverage') || isempty(opts.minCoverage), opts.minCoverage = 0.2; end
if ~isfield(opts,'sizes') || isempty(opts.sizes)
    s = 4:2:nPatients;
    if isempty(s) || s(end) ~= nPatients
        s = [s, nPatients];
    end
    opts.sizes = s;
end
if ~isfield(opts,'nRep') || isempty(opts.nRep),           opts.nRep = 10;         end
if ~isfield(opts,'corrType') || isempty(opts.corrType),   opts.corrType = 'Spearman'; end
if ~isfield(opts,'seed') || isempty(opts.seed),           opts.seed = 42;         end
if ~isfield(opts,'tol') || isempty(opts.tol),             opts.tol = 0.05;        end
end

% ======================================================================
function [statVec, maskIdx, ok] = ea_unified_modelstability_fitone(Xsub, ysub, opts)
% Coverage mask recomputed within THIS subsample (gotcha: never reuse a
% globally-computed mask), then the dense per-element statistic.
validCount = sum(Xsub ~= 0 & ~isnan(Xsub), 2);
coverage   = validCount / size(Xsub,2);
maskIdx    = find(coverage >= opts.minCoverage);

ok = true;
statVec = [];

if isempty(maskIdx)
    ok = false; return
end
if all(isnan(ysub)) || std(ysub,'omitnan') == 0
    ok = false; return
end

Xm = Xsub(maskIdx,:);

if isempty(opts.statFcn)
    statVec = ea_unified_modelstability_defaultcorr(Xm, ysub, opts.corrType);
else
    statVec = opts.statFcn(Xm, ysub);
    statVec = statVec(:);
end

if all(isnan(statVec))
    ok = false;
end
end

% ======================================================================
function r = ea_unified_modelstability_defaultcorr(Xsub, ysub, corrType)
% Vectorized per-row correlation of Xsub (V x n) against ysub (n x 1): a
% single centre-and-normalise pass plus one matrix multiply, never a loop
% over the V (potentially millions of) elements.
%
% For 'Spearman', ranks are ordinal (ties broken by original order via a
% double-sort trick), not the tie-averaged ranks MATLAB's own tiedrank/corr
% would give -- an exact vectorized tie-averaged rank across millions of
% rows at once is expensive, and ties are rare for continuous overlap/
% E-field statistics, so this is an accepted approximation for a stability
% diagnostic (not a p-value engine). Use a custom opts.statFcn wrapping
% corr(...,'Rows','pairwise') if exact tie handling matters to you.
if strcmpi(corrType, 'Spearman')
    [~, ord]  = sort(Xsub, 2);
    [~, M]    = sort(ord, 2);
    [~, ordy] = sort(ysub);
    [~, yv]   = sort(ordy);
else
    M  = Xsub;
    yv = ysub;
end

Mc = M - mean(M, 2);
yc = yv(:)' - mean(yv);

denomY = sqrt(sum(yc.^2));
if denomY == 0
    r = nan(size(Xsub,1), 1);
    return
end

num    = Mc * yc';
denomX = sqrt(sum(Mc.^2, 2));

r = num ./ (denomX * denomY);
r(denomX == 0) = nan;
end

% ======================================================================
function [selPos, selNeg, selPool] = ea_unified_modelstability_select(statVec, maskIdx, opts)
% Selection rule applied UNIFORMLY regardless of how the caller's main
% analysis thresholds its display -- needed for Dice/centroid to be
% comparable across models with different underlying set sizes.
valid = ~isnan(statVec);
vals  = statVec(valid);
ids   = maskIdx(valid);

selPos = []; selNeg = [];

if any(strcmp(opts.tails, {'positive','both'}))
    posMask = vals > 0;
    selPos = ea_unified_modelstability_selecttail(ids(posMask), vals(posMask), opts, true);
end
if any(strcmp(opts.tails, {'negative','both'}))
    negMask = vals < 0;
    selNeg = ea_unified_modelstability_selecttail(ids(negMask), vals(negMask), opts, false);
end

selPool = union(selPos, selNeg);
end

function sel = ea_unified_modelstability_selecttail(ids, vals, opts, isPos)
if isempty(ids)
    sel = [];
    return
end
if isPos
    [~, ord] = sort(vals, 'descend');
else
    [~, ord] = sort(vals, 'ascend'); % most negative first
end
idsSorted = ids(ord);

switch lower(opts.selection)
    case 'topk'
        n = min(opts.k, numel(idsSorted));
    case 'percentile'
        n = max(1, round(opts.percentile/100 * numel(idsSorted)));
    otherwise
        ea_error('ea_unified_modelstability: unknown opts.selection "%s" (expected ''topk'' or ''percentile'').', opts.selection);
end
sel = idsSorted(1:n);
end

% ======================================================================
function r = ea_unified_modelstability_spearman(statA, idxA, statB, idxB)
% Dense-statistic-vector agreement over the elements covered by BOTH
% models (their coverage masks can differ between subsamples).
[~, ia, ib] = intersect(idxA, idxB);
if numel(ia) < 3
    r = nan; return
end
a = statA(ia); b = statB(ib);
valid = ~isnan(a) & ~isnan(b);
if sum(valid) < 3
    r = nan; return
end
r = corr(a(valid), b(valid), 'Type', 'Spearman');
end

% ======================================================================
function [dPos, dNeg, dPool] = ea_unified_modelstability_dice(prevPos, prevNeg, prevPool, curPos, curNeg, curPool)
dPos  = ea_unified_modelstability_dicecoef(prevPos,  curPos);
dNeg  = ea_unified_modelstability_dicecoef(prevNeg,  curNeg);
dPool = ea_unified_modelstability_dicecoef(prevPool, curPool);
end

function d = ea_unified_modelstability_dicecoef(A, B)
% Selected sets are original element indices (row IDs into X), which stay
% identity-stable across subsamples, so intersect/union need no re-
% alignment the way the Spearman comparison above does.
if isempty(A) && isempty(B)
    d = nan; return
end
inter = numel(intersect(A, B));
d = 2*inter / (numel(A) + numel(B));
end

% ======================================================================
function d = ea_unified_modelstability_centroiddist(coords, prevPool, curPool)
if isempty(prevPool) || isempty(curPool)
    d = nan; return
end
c1 = mean(coords(prevPool,:), 1);
c2 = mean(coords(curPool,:), 1);
d = norm(c1 - c2);
end

% ======================================================================
function pt = ea_unified_modelstability_stabpoint(sizes, curve, tol)
% Earliest sample size after which the curve never again leaves a +/-tol
% band around its own final value. No exponential fitting, no elbow
% detection -- exactly the criterion used by the source paper.
valid = find(~isnan(curve));
if numel(valid) < 2
    pt = nan; return
end

finalVal = curve(valid(end));
band = sort([finalVal*(1-tol), finalVal*(1+tol)]);
lo = band(1); hi = band(2);
if lo == hi % finalVal == 0 exactly -- fall back to an absolute band
    lo = -tol; hi = tol;
end

lastViol = 0;
for k = 1:numel(valid)-1 % the final point trivially satisfies its own band
    v = curve(valid(k));
    if v < lo || v > hi
        lastViol = valid(k);
    end
end

if lastViol == 0
    pt = sizes(valid(1));
else
    idxAfter = valid(find(valid > lastViol, 1, 'first'));
    pt = sizes(idxAfter);
end
end

% ======================================================================
function loo = ea_unified_modelstability_loo(X, y, opts)
% Free LOO tier: n leave-one-out models, reused for (a) pairwise agreement
% among all of them (stability at n-1, no extra fitting) and (b) a per-
% patient influence score.
n = size(X,2);

[fullStat, fullMaskIdx, fullOk] = ea_unified_modelstability_fitone(X, y, opts);

stats = cell(n,1);
masks = cell(n,1);
okv   = false(n,1);

for j = 1:n
    keep = true(1,n); keep(j) = false;
    [s, m, ok] = ea_unified_modelstability_fitone(X(:,keep), y(keep), opts);
    stats{j} = s; masks{j} = m; okv(j) = ok;
end

selPos = cell(n,1); selNeg = cell(n,1); selPool = cell(n,1);
for j = 1:n
    if okv(j)
        [selPos{j}, selNeg{j}, selPool{j}] = ea_unified_modelstability_select(stats{j}, masks{j}, opts);
    end
end

haveCoords = ~isempty(opts.coords);
pairR = nan(n,n); pairD = nan(n,n); pairC = nan(n,n);

for a = 1:n-1
    if ~okv(a), continue; end
    for b = a+1:n
        if ~okv(b), continue; end
        pairR(a,b) = ea_unified_modelstability_spearman(stats{a}, masks{a}, stats{b}, masks{b});
        [~,~,pairD(a,b)] = ea_unified_modelstability_dice(selPos{a}, selNeg{a}, selPool{a}, selPos{b}, selNeg{b}, selPool{b});
        if haveCoords
            pairC(a,b) = ea_unified_modelstability_centroiddist(opts.coords, selPool{a}, selPool{b});
        end
    end
end

loo = struct;
loo.pairwise.spearman_mean = mean(pairR(:), 'omitnan');
loo.pairwise.spearman_sd   = std(pairR(:), 'omitnan');
loo.pairwise.dice_mean     = mean(pairD(:), 'omitnan');
loo.pairwise.dice_sd       = std(pairD(:), 'omitnan');
if haveCoords
    loo.pairwise.centroid_mean = mean(pairC(:), 'omitnan');
    loo.pairwise.centroid_sd   = std(pairC(:), 'omitnan');
else
    loo.pairwise.centroid_mean = nan;
    loo.pairwise.centroid_sd   = nan;
end

score = nan(n,1);
if fullOk
    for j = 1:n
        if okv(j)
            score(j) = 1 - ea_unified_modelstability_spearman(fullStat, fullMaskIdx, stats{j}, masks{j});
        end
    end
end

allIdx    = (1:n)';
validMask = ~isnan(score);
[sortedValid, ordValid] = sort(score(validMask), 'descend');
validIdxList   = allIdx(validMask);
invalidIdxList = allIdx(~validMask);

loo.influence.patientIdx = [validIdxList(ordValid); invalidIdxList];
loo.influence.score      = [sortedValid; nan(numel(invalidIdxList),1)];
end

% ======================================================================
function v = ea_unified_modelstability_firstout(fcn)
% Convenience helper for README-style statFcn wrappers -- picks the first
% output of a multi-output call. Not used internally.
[v] = fcn();
end
