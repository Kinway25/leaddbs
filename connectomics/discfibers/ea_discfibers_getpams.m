function pamlist = ea_discfibers_getpams(obj)
% Return list of VATs

% For multiple protocols, we should look for indexed files in the stim
% folder, but mark that they are from the same patient
numPatient = length(obj.allpatients);
pamlist = cell(numPatient,2);   % no mirroring

% here we will add the missing ones (you just need to know how much you have in total, iterate, set to zero for no match)

disp('Construct PAM list...')

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,vta_name,~] = fileparts(obj.M.patient.list{sub});
    parts = strsplit(vta_name,'_');
    subj_tag = parts{1};
    if contains(vta_name,'right')
        AmpTrajDepth = parts{end};
    else
        AmpTrajDepth = parts{end-1};
    end
    subSimPrefix = ['sub-',subj_tag, '_sim-'];
    fiberActivation_file_rh = [subSimPrefix,'fiberActivation_model-ossdbs_hemi-R'];
    fiberActivation_file_lh = [subSimPrefix,'fiberActivation_model-ossdbs_hemi-L'];
    
    leadset_path = '/media/netstim/Konstantin/JRBIDS/derivatives/leaddbs/';

    fiberActivation_merged_rh = [leadset_path,'sub-',subj_tag, '/stimulations/native/gs_20230213012911/Results_rh_sorted_biphasic',filesep,AmpTrajDepth,filesep,fiberActivation_file_rh,'.mat'];
    fiberActivation_merged_lh = [leadset_path,'sub-',subj_tag, '/stimulations/native/gs_20230213012911/Results_lh_sorted_biphasic',filesep,AmpTrajDepth,filesep,fiberActivation_file_lh,'.mat'];

    % only assign one side
    if contains(vta_name,'right')
        pamlist{sub,1} = fiberActivation_merged_rh;
        pamlist{sub,2} = 'skip';
    else
        pamlist{sub,2} = 'skip';
        pamlist{sub,1} = fiberActivation_merged_lh;
    end
end

end
