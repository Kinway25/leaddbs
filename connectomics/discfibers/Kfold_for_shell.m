function [training_sets,test_sets] = Kfold_for_shell(obj,patientsel,threshold_STN_bin)

PT_names = {'MNI/1cJTEbdxV6EKuefRmS6bQ6a'
'MNI/22b959597732'
'MNI/2kZBXs9Tre3o9NAiQWQCb8'
'MNI/4d1f5973eb20'
'MNI/4eVfR2fjrSw2tMQLf9BbQt'
'MNI/8aPF46swc33pKnC6Rs4RaS'
'MNI/8c9c1101b527'
'MNI/91VmEGp8Z8x1ebK6N8kE1h'
'MNI/9gq3BdcisR2T2qGvtjhf6Z'
'MNI/SBN0BBDF'
'MNI/SBN0DS64'
'MNI/SBN1S74E'
'MNI/SBN39922'
'MNI/SBN4E2BE'
'MNI/SBN6F007'
'MNI/SBN9974A'
'MNI/SBN9F7S6'
'MNI/SBNA396B'
'MNI/SBNA94FS'
'MNI/SBNB0SBB'
'MNI/SBNB304A'
'MNI/SBNB40B6'
'MNI/SBNC46CB'
'MNI/SBNCB20B'
'MNI/SBND0360'
'MNI/SBNE1A7F'
'MNI/SBNEA29F'
'MNI/SBNEBCA3'
'MNI/SBNFSSFA'
'MNI/SBNS2366'
'MNI/b39f3f6874d3'
'MNI/beaucMTKeenmkGCfPCx4AQ'
'MNI/c9LM7G63AsUvxvv1T4cj8z'
'MNI/cFfPeBrdMafcN4TFCUExGE'
'MNI/eLQ1iuxG8nQDBr9yreUiGq'
'MNI/f1317979a99e'
'MNI/kY9M6btomZ9jELQQ8WA1wq'
'MNI/mvwcN7XHeFRkCkMeNVuupT'
'MNI/oKcxih6DX58fDpzHxzN1y3'
'MNI/ogqgsgAAsgKi2NAv3iqegW'
'MNI/sMKbfUxmYtXmA2q6KxRveR'
'MNI/uXSVxdRjFN5Mkmz3jbxyWV'};


%N_PTs = length(PT_names);
my_indices = randperm(length(PT_names));
N_folds = 5;

% here you need to split into nearly equal folds, e.g
%N_elements = [4,4,4,4,4,4,4,4,5,5];
N_elements = [8,8,8,9,9];

fold_out = cell(1, N_folds);

gl_counter = 1;

% initiate training and test in the subcohort space
training_sets = logical(zeros(length(threshold_STN_bin),N_folds));
test_sets = logical(zeros(length(threshold_STN_bin),N_folds));

for fold_i = 1:N_folds
    fold_out{1,fold_i} = cell(N_elements(fold_i),1);
    for pt = 1:N_elements(fold_i)
        fold_out{1,fold_i}{pt} = PT_names{my_indices(pt+gl_counter-1)};
        %
        disp(fold_out{1,fold_i}{pt})
    end


    % now you iterate over all VTAs and unselect that patient
    % in both training and test
    for vta_j = 1:length(patientsel)
        temp2=strsplit(obj.M.patient.list{patientsel(vta_j)},'_');
        pt_ID2 = temp2{2}; 

        % check if patient is in this fold
        if ~any(strcmp(fold_out{1,fold_i},pt_ID2))

            %continue
            % for training we need additionally check if it is
            % threshold STN value 
            if threshold_STN_bin(1,vta_j) == 1
                training_sets(vta_j,fold_i) = 1;
            end
        else
            % if yes, select for test
            test_sets(vta_j,fold_i) = 1;
        end
    end


    gl_counter = gl_counter + N_elements(fold_i);
    disp(sum(training_sets(:,fold_i)))
    disp(sum(test_sets(:,fold_i)))
    disp("_______")
end
