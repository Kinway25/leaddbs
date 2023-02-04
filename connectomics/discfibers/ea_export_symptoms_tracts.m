function ea_export_symptoms_tracts(obj, vals_threshold, negative_vals)

% this function creates symptom-specific tracts from FF
% Important: use either positive or negative tracts
% Negative tracts can be used for soft sode-effects (their weights will be abs(norm((vals)))


% we want to define the model outside of cross-validation
if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patientsel = obj.patientselection;
end

% get fiber model (vals) and corresponding indices (usedidx)
if obj.cvlivevisualize
    [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel);
    obj.draw(vals,fibcell,usedidx)
    drawnow;
else
    [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel);
end

% vals_threshold - only fibers with vals above will be taken (atm, hardcoded to 0.25)
vals_threshold = 0.25;

% if negative tracts, make a note in the pathway name
% those are for soft-threshold side-effects
negative_vals = 0;

% unlike export fibscore model, we do not need to map to global space,
% but only have the global indices



for voter=1:size(vals,1)  % I would restrict to one voter for now
    for side = 1:size(vals,2)

        % vals_threshold has to be used here

        gl_indices = obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{voter,side});

        % normalized to 0-1 if necessary (apply abs if negative)
        if min(vals{voter,side}) < 0 && max(vals{voter,side}) > 0
            disp("Choose only positive or negative tracts!")
            return;
        elseif min(vals{voter,side}) < 0
            vals{voter,side} = abs(vals{voter,side});
        end

        if max(vals{voter,side}) > 1.0
            vals{voter,side} = vals{voter,side} / max(vals{voter,side});
        end


        % Let's reverse the order:
        % 1) Group based on Corr. matrix 
        % 2) Then bin using Jenks Natural Breaks, more bins for higher vals

        % load cfile
        [filepath,~,~] = fileparts(obj.leadgroup);
        if obj.multi_pathways == 1
            cfile = [filepath,filesep,'merged_pathways.mat'];
        else
            cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
        end
        load(cfile, 'fibers', 'idx');

        % we can check and remove fibers that are too far away from the
        % target (but track the indices!)


        % create fibers and idx only from remaining fibers
        % but remember to use gl_indices
        fibers_with_vals = [];
        fibers_orig_ind = [];
        idx_with_vals = idx(gl_indices);
        fibers_all = fibers;
        idx_all = idx;

        for i = 1:length(gl_indices)
            fibers_with_vals = cat(1,fibers_with_vals,fibers_all(fibers_all(:,4) == gl_indices(i),:));
            fibers_with_vals(fibers_with_vals(:,4) == gl_indices(i),4) = i;
            fibers_orig_ind = cat(1,fibers_orig_ind,fibers_all(fibers_all(:,4) == gl_indices(i),4));
        end


        % Min Jae's magic to get PRX_M
        % tracts should be ordered exactly like in gl_indices
        sig = 1.5; % Gaussian Kernel Size 
        if side == 1
            trgt_coor = [11.7631,-13.9674,-8.85935]; % target coordinate (e.g STN), right and left
        else
            trgt_coor = [-11.7631,-13.9674,-8.85935];
        end
            
        vcnty_thr = 15; % 15 mm works well. threshold vicinity mask | N.B: units are in mm
        spatio_corr_mat = fiber_spatial_rev3(fibers_with_vals,trgt_coor,sig,vcnty_thr);


        % test matrix
%         d = rand(length(gl_indices),1); % The diagonal values
%         t = triu(bsxfun(@min,d,d.').*rand(length(gl_indices)),1); % The upper trianglar random values
%         M = diag(d)+t+t.'; % Put them together in a symmetric matrix
% 
%         PRX_M = triu(M,1);
%         PRX_M = triu(euc_dist_mat,1);
        %PRX_M = triu(spatio_corr_mat(1:100,1:100),1);
        PRX_M = triu(spatio_corr_mat,1)';
        % hierarchical clustering of the correlation matrix
        
        % you can check mean(PRX_M) to decide on the clustering threhsold

        prox_cluster_threshold = 0.1; % should be associated with some meaningful distance metric
                                      % the larger the value, the more groups you get

        %PRX_M = PRX_M - eye(size(PRX_M));

        % and convert to a vector (as pdist)
        dissimilarity = 1 - PRX_M(find(PRX_M))';
        %check_matrix = squareform(dissimilarity);

        % decide on a cutoff
        % remember that 0.4 corresponds to corr of 0.6!
        cutoff = 1 - prox_cluster_threshold; 
        
        %% This is wrong: linkage is used on spatially distributed points, not the spatial difference
        % perform complete linkage clustering
        %Z = linkage(dissimilarity,'complete','correlation');
        Z = linkage(dissimilarity);
% 
%         pdist_vector_size = ((size(PRX_M,1) ^ 2)- size(PRX_M,1))/2;
%         PRX_M_pdist = zeros(1,pdist_vector_size); % it is actually quadratic by definition
%         gl_counter = 1;
%         for i = 1:size(PRX_M,1)
%             for j = 1:size(PRX_M,2)
%                 if j < i
%                     PRX_M_pdist(1,gl_counter) = PRX_M(j,i);
%                     gl_counter = gl_counter + 1;
%                 end
%             end
%         end
% 
%         %Z = linkage(PRX_M,'complete','euclidean');
%         Z = linkage(PRX_M_pdist);
        
        % group the data into clusters
        % (cutoff is at a correlation of 0.5)
        groups = cluster(Z,'cutoff',cutoff,'criterion','distance');


        %% we can visualize groups using OSS-DBS output format
        for i = 1:length(idx_with_vals)
            fibers_with_vals(fibers_with_vals(:,4) == i,5) = groups(i)-4;
        end
        origNum = length(idx);
        connectome_name = 'placeholder';

        % rename just to be able to load later
        fibers = fibers_with_vals; % IMPORTANT: remember that 4th column are local indices here!
        idx = idx_with_vals;
        % then save to mat origNum, fibres, connectome_name and idx

        % gl_indices(groups == 1) returns fibers from the first spatially
        % defined pathway

        % usedidx{voter,side}(groups == 1) can be used to retrieve vals?

        if groups > 10
            disp("Warning: number of spatially distributed pathways above 10")
            disp("Number of val bins per pathway will be reduced to 5")
            max_N_pathways = 5; % per spatially defined pathway
        else
            max_N_pathways = 10;
        end


        % Now bin within blocks based on variance 
        pathway_index_list = cell(1,length(unique(groups)));
        pathway_mean_vals = cell(1,length(unique(groups)));

        disp("Number of spatial groups "+string(length(unique(groups))))

        for group_i = 1:length(unique(groups))

            %clc; clear output sub_array;
    
            if length(vals{voter,side}(groups == group_i)) > 1
    
                % check performance with real data
                if max_N_pathways > length(find(groups == group_i))
                    max_N_pathways_bin = length(find(groups == group_i));
                else
                    max_N_pathways_bin = max_N_pathways;
                end
    
                metric_vals = zeros(max_N_pathways_bin-1,1);
                for class_number = max_N_pathways_bin:-1:2
                    
                    input = vals{voter,side}(groups == group_i);  % test
                    idx_group = gl_indices(groups == group_i);
    
                    % Otsu's variance based method
                    [thresh, metric] = multithresh(input,class_number);
                    if metric == 0 && length(input) > 1
                        % did not work, just split by half
                        thresh = mean(input);
                        break
                    end
                    metric_vals(class_number-1) = metric;
                    if class_number == max_N_pathways_bin
                        metric10 = metric;
                    else
                        if metric < 0.75 * metric10 % check this condition empirically, allow 25% reduction
                            fprintf("N classes: %d \n", class_number + 1 )
                            % recompute thresholds
                            [thresh, metric] = multithresh(input,class_number + 1);
                            break
                        end
                    end
                end
    
                %plot(metric_vals,1:length(metric_vals))
    
                pathway_index_list{1,group_i} = cell(1,length(thresh) + 1);
                pathway_mean_vals{1,group_i} = zeros(1,length(thresh) + 1);
    
                for thresh_i = 1:length(thresh) + 1
    
                    if thresh_i == 1
    
                        % e.g. drop pathways with vals < 0.5 and N fibers < 10
                        pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input < thresh(thresh_i)));
    
                        if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(find(idx_group(input < thresh(thresh_i)))) < 10
                            pathway_index_list{1,group_i}{thresh_i} = 0;
                        else
                            pathway_index_list{1,group_i}{thresh_i} = idx_group(input < thresh(thresh_i));
                        end
    
                    elseif thresh_i == length(thresh) + 1
    
                        pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input >= thresh(thresh_i-1)));
                        pathway_index_list{1,group_i}{thresh_i} = idx_group(input >= thresh(thresh_i-1));
    
    %                     if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(pathway_index_list{1,group_i}{thresh_i}) < 10
    %                         pathway_index_list{1,group_i}{thresh_i} = 0;
    %                     else
    %                         pathway_index_list{1,group_i}{thresh_i} = gl_indices(find(input > thresh(thresh_i-1)));
    %                     end
                    else 
                        pathway_mean_vals{1,group_i}(thresh_i) = mean(input(input >= thresh(thresh_i-1) & input < thresh(thresh_i)));
                        
                        % e.g. drop pathways with vals < 0.5 and N fibers < 10
                        if pathway_mean_vals{1,group_i}(thresh_i) < 0.5 && length(find(input >= thresh(thresh_i-1) & input < thresh(thresh_i))) < 10
                            pathway_index_list{1,group_i}{thresh_i} = 0;
                        else
                            pathway_index_list{1,group_i}{thresh_i} = idx_group(input >= thresh(thresh_i-1) & input < thresh(thresh_i));
                        end
    
                    end
                end
            else
                pathway_mean_vals{1,group_i}(1) = vals{voter,side}(groups == group_i);
                pathway_index_list{1,group_i}(1) = gl_indices(groups == group_i);
            end

            % save as a pathway
            %fibers, idx, fourindex, ea_fibformat and fibers_glob_index
            ea_fibformat = '1.0';
            fourindex = 1;

            for i= 1:length(pathway_index_list{1,group_i})
                fibers_pathway = [];
                fibers_glob_ind = [];

                if length(pathway_index_list{1,group_i}) == 1
                    fibers_pathway = fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i},:);
                    fibers_pathway(:,4) = 1;
                    fibers_glob_ind = pathway_index_list{1,group_i};
                    idx_pathway(1) = size(fibers_pathway,1);
                else
                    idx_pathway = zeros(length(pathway_index_list{1,group_i}{1,i}),1);
                    for j = 1:length(pathway_index_list{1,group_i}{1,i})
                        fibers_pathway = cat(1,fibers_pathway,fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),:));
                        fibers_pathway(fibers_pathway(:,4) == pathway_index_list{1,group_i}{1,i}(j),4) = j;
                        fibers_glob_ind = cat(1,fibers_glob_ind,fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),4));
                        idx_pathway(j) = size(fibers_all(fibers_all(:,4) == pathway_index_list{1,group_i}{1,i}(j),4),1);
                    end
                end
                ftr.fibers = fibers_pathway;
                ftr.idx = idx_pathway;
                ftr.ea_fibformat = ea_fibformat;
                ftr.fourindex = fourindex;

                filename = char(compose('pathway_side_%d_%d_val.mat',side,pathway_mean_vals{1,group_i}(i)));
                % you should also add the symptom name
                save([filepath,filesep,filename],'-struct','ftr')
            end

        end

    end
end        










end