function [emp_p, sig] = ea_unified_permutation_threshold(realvals, nullvals, obj, mode)
% Turns a permutation null distribution into empirical p-values and a
% significance mask, for the unified mapping explorer's
% 'Permutation Threshold (...)' multcompstrategy branch. The null
% distribution itself is built separately, by ea_unified_permutation_nulldist.m
% (the parallelizable, expensive part).
%
% realvals - V x 1 vector of observed (signed) test statistics (e.g.
%            r-values), already restricted to the rows that were actually
%            tested.
% nullvals - the null distribution to compare realvals against. Shape/type
%            depends on mode:
%              'uncorrected' - V x Nperm matrix, one null row per feature
%                               (this cell's own ea_unified_permutation_nulldist
%                               output). Each feature is compared one-tailed
%                               against only the tail of its own null that
%                               matches its own sign (positive real value ->
%                               count null >= real value; negative real
%                               value -> count null <= real value), then
%                               doubled for two-tailedness. Comparing
%                               against only the matching tail (rather than
%                               pooling both tails via abs()) means an
%                               asymmetric null -- e.g. this fiber's shuffled
%                               r happening to swing further positive than
%                               negative under skewed outcome data -- can't
%                               make a genuinely extreme value in the other
%                               direction look unremarkable. This is NOT
%                               corrected for testing many features at once
%                               (comparable to the 'Uncorrected' strategy,
%                               just built empirically instead of
%                               parametrically).
%              'maxstat'     - struct with fields .max and .min, each
%                               1 x Nperm, shared by every feature across
%                               all cells/sides/groups (built by
%                               ea_unified_corrsignan.m). .max(p) is the
%                               largest *signed* statistic seen anywhere in
%                               the whole analysis under permutation p;
%                               .min(p) is the smallest (most negative).
%                               A fiber with a positive real value is
%                               tested one-tailed against .max; a fiber
%                               with a negative real value is tested
%                               one-tailed against .min -- each fiber only
%                               ever competes against the tail matching its
%                               own sign, so an asymmetric null (positive
%                               and negative spurious correlations behaving
%                               differently, e.g. under skewed outcome
%                               data) doesn't let one tail contaminate the
%                               other's significance assessment the way
%                               pooling via abs() would. This is the
%                               Nichols & Holmes (2002) signed max-statistic
%                               method, and directly controls the
%                               family-wise error rate (FWER) -- no further
%                               FDR/Bonferroni step needed on top of it.
% obj      - explorer object; uses obj.alphalevel.
% mode     - 'uncorrected' (default) or 'maxstat'.
%
% Returns:
% emp_p - V x 1 empirical two-tailed p-values (NaN where realvals is NaN).
%         For 'maxstat' this is the one-tailed p from whichever of Tmax/Tmin
%         matched the feature's sign (see sig below for why); for
%         'uncorrected' it's the sign-matched one-tailed p already doubled.
% sig   - V x 1 logical, true where significant at obj.alphalevel. For
%         'maxstat' this compares each one-tailed emp_p against
%         obj.alphalevel/2, since two one-tailed tests (one per direction)
%         are being run to keep the overall two-tailed FWER at
%         obj.alphalevel -- the same reasoning as reading a two-tailed
%         t-test off a table using alpha/2 in each tail. For 'uncorrected',
%         emp_p is already the doubled two-tailed p, so it's compared
%         directly against the full obj.alphalevel.

if nargin<4 || isempty(mode)
    mode='uncorrected';
end

switch lower(mode)
    case 'maxstat'
        Nperm = numel(nullvals.max);
        emp_p = ones(size(realvals)); % real value exactly 0 (no evidence either way) stays non-significant
        pos = realvals > 0;
        neg = realvals < 0;
        exceedPos = sum(nullvals.max >= realvals(pos), 2, 'omitnan');
        emp_p(pos) = (exceedPos + 1) ./ (Nperm + 1);
        exceedNeg = sum(nullvals.min <= realvals(neg), 2, 'omitnan');
        emp_p(neg) = (exceedNeg + 1) ./ (Nperm + 1);

        emp_p(isnan(realvals)) = nan;
        sig = emp_p <= (obj.alphalevel/2);
    otherwise
        Nperm = size(nullvals,2);
        exceed = zeros(size(realvals));
        pos = realvals > 0;
        neg = realvals < 0;
        exceed(pos) = sum(nullvals(pos,:) >= realvals(pos), 2, 'omitnan');
        exceed(neg) = sum(nullvals(neg,:) <= realvals(neg), 2, 'omitnan');

        emp_p = min(2 * (exceed + 1) ./ (Nperm + 1), 1); % one-tailed, matching sign, doubled for two-tailedness
        emp_p(~pos & ~neg) = 1; % real value exactly 0 -- no evidence either way

        emp_p(isnan(realvals)) = nan;
        sig = emp_p <= obj.alphalevel;
end

sig(isnan(realvals)) = false;
