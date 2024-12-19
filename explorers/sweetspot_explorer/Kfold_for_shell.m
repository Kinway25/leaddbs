function [training_sets,test_sets] = Kfold_for_shell(obj,patientsel,patientsel_train,threshold_STN_bin)

load('/home/forel/Documents/data/JB_project/JK_SW_table.mat')
PT_names = string(unique(data_flat_SW.subject));
%load('/home/forel/Documents/data/JB_project/JB_SW_table.mat')
%PT_names = unique(data_flat_JB.pt_label);

%N_PTs = length(PT_names);
my_indices = randperm(length(PT_names));
N_folds = 16;

% here you need to split into nearly equal folds, e.g
%N_elements = [4,4,4,4,4,4,4,4,5,5];
N_elements = ones(1,size(PT_names,1));
%N_elements = [8,8,8,9,9];

fold_out = cell(1, N_folds);

gl_counter = 1;

% initiate training and test in the subcohort space
%training_sets = logical(zeros(length(threshold_STN_bin),N_folds));
training_sets = logical(zeros(length(patientsel),N_folds));
%test_sets = logical(zeros(length(threshold_STN_bin),N_folds));
test_sets = logical(zeros(length(patientsel),N_folds));
patientsel_full = [];

for fold_i = 1:N_folds
    fold_out{1,fold_i} = cell(N_elements(fold_i),1);
    for pt = 1:N_elements(fold_i)
        fold_out{1,fold_i}{pt} = PT_names{my_indices(pt+gl_counter-1)};
        disp(fold_out{1,fold_i}{pt})
    end


%     % now you iterate over all VTAs and unselect that patient
%     % in both training and test
%     for vta_j = 1:length(patientsel)
%         temp2=strsplit(obj.M.patient.list{patientsel(vta_j)},'_');
%         pt_ID2 = temp2{2}; 
% 
%         % check if patient is in this fold
%         if ~any(strcmp(fold_out{1,fold_i},pt_ID2))
% 
%             %continue
%             % for training we need additionally check if it is
%             % threshold STN value 
%             if threshold_STN_bin(1,vta_j) == 1
%                 training_sets(vta_j,fold_i) = 1;
%             end
%         end
% %         else
% %             % if yes, select for test
% %             test_sets(vta_j,fold_i) = 1;
% %         end
%     end

    for vta_j = 1:length(patientsel)
        temp2=strsplit(obj.M.patient.list{patientsel(vta_j)},'-');
        pt_ID2 = ['0',temp2{2}(1:2)]; 

        if ismember(vta_j,patientsel_train) && ~any(strcmp(fold_out{1,fold_i},pt_ID2))
            if threshold_STN_bin(1,vta_j) == 1
                training_sets(vta_j,fold_i) = 1;
            end
            patientsel_full = [patientsel_full, vta_j];
        elseif any(strcmp(fold_out{1,fold_i},pt_ID2))
            % if yes, select for test
            test_sets(vta_j,fold_i) = 1;
            patientsel_full = [patientsel_full, vta_j];
        end
    end

    gl_counter = gl_counter + N_elements(fold_i);
    %disp(sum(training_sets(:,fold_i)))
    %disp(sum(test_sets(:,fold_i)))
    %disp("_______")
end