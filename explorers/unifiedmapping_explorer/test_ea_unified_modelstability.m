function test_ea_unified_modelstability()
% test_ea_unified_modelstability  Regression tests for ea_unified_modelstability.m. No
% real Lead-DBS data required -- everything is synthetic. Run with:
%   test_ea_unified_modelstability
% Prints a pass/fail summary and errors if anything failed. Whole suite
% runs in well under 30s.

ticAll = tic;
results = {};

fprintf('Running ea_unified_modelstability test suite...\n\n');

% --- Tests 1+2: null case vs strong-signal case ------------------------
% Both share the same subsample schedule/patient count/seed so the only
% thing that differs is whether X actually contains a subset of elements
% that truly correlates with y. Asserting the RELATIVE ordering (signal
% clearly more stable/selective than null) is more robust than hard-coding
% an absolute threshold for the null case alone, since even under pure
% noise the additive/nested subsample design gives adjacent models
% overlapping patients -- so some baseline agreement is expected by
% construction, not just under a true effect. Dice on the selected set is
% the sharper discriminator here (dense Spearman agreement can rise with n
% under null too, simply because the correlation *estimate*'s own sampling
% noise shrinks as n grows).
nTrue = 40;
[nullS, sigS] = run_null_and_signal(nTrue);

% Loose absolute sanity bound only -- empirically, even pure null noise
% reaches ~0.6 pooled Dice by n=24 here (adjacent additive-subsample
% models share most of their patients, so top-k membership churns less
% than independent resampling would suggest). The relative check below is
% the real discriminator; this just guards against Dice saturating near 1
% with no true structure at all, which WOULD indicate a bug.
results{end+1} = check( ...
    'null case does not saturate near perfect agreement (Dice)', ...
    isfinite(nullS.dice_mean.pooled(end)) && nullS.dice_mean.pooled(end) < 0.75, ...
    sprintf('final pooled Dice = %.3f (expected < 0.75)', nullS.dice_mean.pooled(end)));

results{end+1} = check( ...
    'signal case is clearly more stable than null (Dice)', ...
    (sigS.dice_mean.pooled(end) - nullS.dice_mean.pooled(end)) > 0.2, ...
    sprintf('signal Dice=%.3f, null Dice=%.3f', sigS.dice_mean.pooled(end), nullS.dice_mean.pooled(end)));

results{end+1} = check( ...
    'signal case stability point is finite', ...
    isfinite(sigS.stabilityPoint.spearman), ...
    sprintf('stabilityPoint.spearman = %s', mat2str(sigS.stabilityPoint.spearman)));

recall = numel(intersect(sigS.finalModel.selPos, 1:nTrue)) / nTrue;
results{end+1} = check( ...
    'true elements dominate the final selection', ...
    recall > 0.7, ...
    sprintf('recall of true positive elements in final top-k = %.2f (expected > 0.7)', recall));

% --- Test 3: coverage-mask regression -----------------------------------
% Pad the null case with a large block of all-zero elements. If the
% within-subsample coverage mask is doing its job, those elements are
% never covered (0 patients with a non-zero value in any subsample) and
% get excluded before the statistic is ever computed on them -- so the
% Spearman curve should be UNCHANGED (not just similar) versus the
% unpadded null case, since the patient subsample ordering (opts.seed) and
% the original 2000 elements are generated identically in both runs.
nullPadded = run_null_padded();
maxDiff = max(abs(nullS.spearman_mean - nullPadded.spearman_mean), [], 'omitnan');
results{end+1} = check( ...
    'coverage mask excludes all-zero padding (Test 1 vs padded)', ...
    maxDiff < 1e-8, ...
    sprintf('max abs diff in spearman_mean = %.2e (expected < 1e-8)', maxDiff));

% --- Test 4: both modalities (with/without coords) ----------------------
[Swith, Swithout] = run_both_modalities();
results{end+1} = check( ...
    'fiber-like run (coords supplied) produces a finite centroid metric', ...
    ~isempty(Swith.centroid_mean) && any(isfinite(Swith.centroid_mean)), ...
    'S.centroid_mean was empty or all-NaN');

results{end+1} = check( ...
    'voxel-like run (no coords) skips the centroid metric cleanly', ...
    isempty(Swithout.centroid) && isempty(Swithout.centroid_mean) && isnan(Swithout.stabilityPoint.centroid), ...
    'S.centroid/S.centroid_mean were not empty, or stabilityPoint.centroid was not NaN');

% --- Summary -------------------------------------------------------------
fprintf('\n--- Summary ---\n');
nPass = 0;
for i = 1:numel(results)
    r = results{i};
    if r.pass
        nPass = nPass + 1;
        fprintf('[PASS] %s\n', r.name);
    else
        fprintf('[FAIL] %s -- %s\n', r.name, r.message);
    end
end
fprintf('%d/%d tests passed. Elapsed: %.1fs\n', nPass, numel(results), toc(ticAll));

if nPass < numel(results)
    error('test_ea_unified_modelstability:failed', '%d test(s) failed.', numel(results)-nPass);
end
end

% ======================================================================
function r = check(name, pass, message)
r = struct('name', name, 'pass', logical(pass), 'message', message);
end

% ======================================================================
function opts = common_opts(nPatients, k)
opts = struct('sizes', 4:4:nPatients, 'nRep', 5, 'seed', 99, ...
    'k', k, 'tails', 'positive', 'minCoverage', 0.2);
end

% ======================================================================
function [nullS, sigS] = run_null_and_signal(nTrue)
nPatients = 24; nElements = 2000;

rng(1);
yNull = randn(nPatients,1);
Xnull = randn(nElements, nPatients);
nullS = ea_unified_modelstability(Xnull, yNull, common_opts(nPatients, nTrue));

rng(2);
ySig = randn(nPatients,1);
Xsig = randn(nElements, nPatients);
Xsig(1:nTrue,:) = 1.2*ySig' + 1.0*randn(nTrue, nPatients);
sigS = ea_unified_modelstability(Xsig, ySig, common_opts(nPatients, nTrue));
end

% ======================================================================
function nullPadded = run_null_padded()
nPatients = 24; nElements = 2000;

rng(1); % identical seed/generation order as the null case above
yNull = randn(nPatients,1);
Xnull = randn(nElements, nPatients);
Xpadded = [Xnull; zeros(5000, nPatients)];

nullPadded = ea_unified_modelstability(Xpadded, yNull, common_opts(nPatients, 40));
end

% ======================================================================
function [Swith, Swithout] = run_both_modalities()
nPatients = 20; nElements = 500;

rng(3);
y = randn(nPatients,1);
X = randn(nElements, nPatients);
X(1:30,:) = 1.0*y' + 1.0*randn(30, nPatients);

optsWith = common_opts(nPatients, 30);
optsWith.coords = randn(nElements,3)*10; % arbitrary mm-like coordinates
Swith = ea_unified_modelstability(X, y, optsWith);

optsWithout = common_opts(nPatients, 30); % coords omitted
Swithout = ea_unified_modelstability(X, y, optsWithout);
end
