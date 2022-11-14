function permpRank = ea_compute_pval_SeqRankAgr(Ihat_iter)

    %%  Evaluating sequential rank agreement between shuffles

    % Implemented following Ekstrom et al
    % https://doi.org/10.1093/biostatistics/kxy017



    % rank all data columnwise
    rankedScores = zeros(size(Ihat_iter{1,1},1),size(Ihat_iter,2));
    for i = 1:size(rankedScores,2)
        [~,rankedScores(:,i)]  = ismember(Ihat_iter{1,i},unique(Ihat_iter{1,i}));
    end
    cardinality = size(unique(Ihat_iter{1,i}),1); % check if ties in some
    
    % compute agreement
    SeqRankAgr = ea_compute_SeqRankAgr(rankedScores,cardinality);
    
    % permute columnwise
    disp('Computing Sequential Rank Agreement for permuted ranks')
    if length(rankedScores(:,1)) <= 7
        num_perms = size(perms(rankedScores(:,1)),1);
    else
        num_perms = 25000;
    end
    
    % compute agreement for each permutation
    seqRankAgr_perm = zeros(num_perms,1);

%     % first generate permutations for each iteration separately
%     % This part is wrong! We reshuffle vectors separately, but this should
%     % be implemented differently
%     shuffledRanks = cell(size(rankedScores,2),1);
%     for i = 1:size(rankedScores,2)
%         shuffledRanks{i} = ea_shuffle(rankedScores(:,i),num_perms,1:length(rankedScores(:,i)),i*10); % have to provide the rngseed here
%     end

%     % comlete random permutation
%     for p = 1:num_perms
%         % combine permuted iteration ranks for p
% 
%         perm_rankedScores = zeros(size(rankedScores));
%         for i = 1:size(rankedScores,2)
%             %perm_rankedScores(:,i) = shuffledRanks{i,1}(p,:)';
%             perm_rankedScores(:,i) = randperm(size(rankedScores,1))';
%         end
% 
%         seqRankAgr_perm(p) = ea_compute_SeqRankAgr(perm_rankedScores,cardinality);
%     end

   % one element per iteration permutation
    for p = 1:num_perms
        % combine permuted iteration ranks for p

        % permute locally (just one entry)
        perm_rankedScores = rankedScores;
        for i = 1:size(rankedScores,2)
            %perm_rankedScores(:,i) = shuffledRanks{i,1}(p,:)';
            swapidx = randperm(numel(perm_rankedScores(:,i)), 3);
            a = perm_rankedScores(:,i);
            a(swapidx) = a(fliplr(swapidx));
            perm_rankedScores(:,i) = a;
        end

        seqRankAgr_perm(p) = ea_compute_SeqRankAgr(perm_rankedScores,cardinality);
    end


    % check the percentile of the original SeqRankAgr
    Count = sum(seqRankAgr_perm > SeqRankAgr);
    permpRank = 1.0 - Count / length(seqRankAgr_perm);



function SeqRankAgr = ea_compute_SeqRankAgr(rankedScores,card)
     
SampleStdErr = zeros(size(rankedScores,1),1);
SeqRankAgr_nom = 0;
for j = 1:size(rankedScores,1)
    nominator = 0;
    for i = 1:size(rankedScores,2)
        nominator = nominator + (rankedScores(j,i)-mean(rankedScores(j,:))).^2;
    end
    SampleStdErr(j) = sqrt(nominator / (size(rankedScores,2)- 1));

    SeqRankAgr_nom = SeqRankAgr_nom + (size(rankedScores,2)- 1)*SampleStdErr(j).^2;
end
SeqRankAgr = sqrt(SeqRankAgr_nom / ((size(rankedScores,2)- 1) * card));

