function [training_sets,test_sets] = LOPO_Cologne(obj,patientsel)

load('/home/forel/Documents/data/JB_project/Cologne/LOPOtableCologne_2mA.mat')
PT_names = string(unique(data_flat.subject));

% load('/home/forel/Documents/data/JB_project/JK_SW_table_18.mat')
% PT_names = string(unique(data_flat_SW.subject));

%N_PTs = length(PT_names);
my_indices = randperm(length(PT_names));
N_folds = length(PT_names);

N_elements = ones(1,N_folds);
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



    for vta_j = 1:length(patientsel)
        temp2=strsplit(obj.M.patient.list{patientsel(vta_j)},'-');
        pt_ID2 =  ['0',temp2{2}(5:6)]; 

        if ~any(strcmp(fold_out{1,fold_i},pt_ID2))
            training_sets(vta_j,fold_i) = 1;
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