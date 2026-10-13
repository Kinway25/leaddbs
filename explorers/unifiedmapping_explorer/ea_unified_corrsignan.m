function vals=ea_unified_corrsignan(vals,ps,obj,permdata)
% Unified-mapping-specific copy of ea_corrsignan.m. Behaves identically for
% 'FDR'/'Bonferroni'/uncorrected. Adds two permutation-based strategies,
% which need the raw per-fiber data (not just the already-computed
% p-values) to build null distributions -- see
% ea_unified_permutation_nulldist.m / ea_unified_permutation_threshold.m.
% permdata is prepared by ea_unifiedmapping_calcstats.m alongside vals/ps.
%
% obj.multcompstrategy (case-insensitive, matched by substring so exact
% GUI label wording doesn't matter):
%   contains 'permutation' & 'max'  -> 'Permutation Threshold (max-statistics)':
%       signed max-statistic (Nichols & Holmes 2002) permutation test. Two
%       null distributions are shared across every fiber/voxel in every
%       group/side cell: the largest signed statistic seen anywhere per
%       permutation (Tmax), and the smallest/most negative (Tmin). A fiber
%       is tested one-tailed against whichever of the two matches its own
%       sign (at alpha/2, so the combined two-tailed FWER stays at
%       obj.alphalevel) -- this avoids letting an asymmetric null (e.g.
%       from skewed outcome data) contaminate one direction's significance
%       assessment with the other direction's extremes, which pooling via
%       abs() into a single null would do. No further FDR/Bonferroni
%       needed on top of it.
%   contains 'permutation' only     -> 'Permutation Threshold (Uncorr)':
%       each feature gets its own null distribution from its own
%       permutations. This is NOT corrected for testing many features at
%       once (comparable to the 'Uncorrected' strategy, just built
%       empirically instead of parametrically).

if nargin<4
    permdata=[];
end

strategy=lower(obj.multcompstrategy);

if startsWith(strategy,'permutation') && ~isempty(permdata)
    if contains(strategy,'max')
        permmode='maxstat';
    else
        permmode='uncorrected';
    end

    % Build every cell's null distribution first (needed either way, and
    % the max-statistic mode additionally needs all of them collected
    % before it can build the single shared null below).
    nulldist=cell(size(vals));
    for cellentry=1:numel(vals)
        entry=permdata{cellentry};
        if isempty(entry) || isempty(entry.nonemptyidx)
            continue
        end
        nulldist{cellentry}=ea_unified_permutation_nulldist(entry,obj);
    end

    if strcmp(permmode,'maxstat')
        validcells=~cellfun(@isempty,nulldist);
        if any(validcells)
            allnull=cat(1,nulldist{validcells});
            sharednull.max=max(allnull,[],1,'omitnan'); % Tmax: largest signed statistic anywhere, per permutation
            sharednull.min=min(allnull,[],1,'omitnan'); % Tmin: smallest (most negative) signed statistic anywhere, per permutation
        end
    end

    for cellentry=1:numel(vals)
        entry=permdata{cellentry};
        if isempty(entry) || isempty(entry.nonemptyidx)
            continue
        end
        realvals=vals{cellentry}(entry.nonemptyidx);
        if strcmp(permmode,'maxstat')
            [~,sig]=ea_unified_permutation_threshold(realvals,sharednull,obj,'maxstat');
        else
            [~,sig]=ea_unified_permutation_threshold(realvals,nulldist{cellentry},obj,'uncorrected');
        end
        keep=false(size(vals{cellentry}));
        keep(entry.nonemptyidx)=sig;
        vals{cellentry}(~keep)=nan;
    end
    return
end

allvals=cat(1,vals{:});
allps=cat(1,ps{:});

nnanidx=~isnan(allvals);
numtests=sum(nnanidx);

switch strategy
    case 'fdr'
        pnnan=allps(nnanidx); % pvalues from non-nan value entries
        [psort,idx]=sort(pnnan);
        pranks=zeros(length(psort),1);
        for rank=1:length(pranks)
            pranks(idx(rank))=rank;
        end
        pnnan=pnnan.*numtests;
        pnnan=pnnan./pranks;
        allps(nnanidx)=pnnan;
    case 'bonferroni'
        allps(nnanidx)=allps(nnanidx).*numtests;
end

allps(~nnanidx)=1; % set p values from values with nan to 1
allvals(allps>obj.alphalevel)=nan; % delete everything nonsignificant.

% feed back into cell format:
cnt=1;
for cellentry=1:numel(vals)
    vals{cellentry}(:)=allvals(cnt:cnt+length(vals{cellentry})-1);
    ps{cellentry}(:)=allps(cnt:cnt+length(vals{cellentry})-1);
    cnt=cnt+length(vals{cellentry});
end
