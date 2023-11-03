function pamlist = ea_discfibers_getpams(obj)
% Return list of VATs
disp('Construct PAM list...')

PAM_mirror_enabled = 1;

if PAM_mirror_enabled == 1
    ea_warndlg("PAM mirroring is used, make sure the connecome is index-symmetric!")

    numPatient = length(obj.allpatients);
    pamlist = cell(numPatient*2,2);  

else
    % For multiple protocols, we should look for indexed files in the stim
    % folder, but mark that they are from the same patient
    numPatient = length(obj.allpatients);
    pamlist = cell(numPatient,2);   % no mirroring
end
% here we will add the missing ones (you just need to know how much you have in total, iterate, set to zero for no match)

%pamlist = cell(numPatient,2);   % custom force to one side (aka PseudoM)

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    subSimPrefix = [subj_tag, '_sim-'];

    %[~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    %subj_tag = subj_tag(1:end-31);
    %subSimPrefix = [subj_tag, '_sim-'];

%     if sub <= 35
%         pamlist{sub,1} = [obj.allpatients{sub}(1:76),subj_tag,'/stimulations/native/gs_06stimfit/',subSimPrefix,'fiberActivation_model-ossdbs_hemi-R.mat'];
%     elseif sub > 35 && sub <= 70
%         pamlist{sub,1} = [obj.allpatients{sub}(1:76),subj_tag,'/stimulations/native/gs_20230506011354/',subSimPrefix,'fiberActivation_model-ossdbs_hemi-R.mat'];
%     elseif sub > 70 && sub <= 105
%         pamlist{sub,1} = [obj.allpatients{sub}(1:76),subj_tag,'/stimulations/native/gs_06stimfit/',subSimPrefix,'fiberActivation_model-ossdbs_hemi-L.mat'];
%     else 
%         pamlist{sub,1} = [obj.allpatients{sub}(1:76),subj_tag,'/stimulations/native/gs_20230506011354/',subSimPrefix,'fiberActivation_model-ossdbs_hemi-L.mat'];
%     end

    actual_location = ['/media/konstantin/ba/_work/BIDSdataPAM/derivatives/leaddbs/',obj.allpatients{sub}(48:end)];
    pamlist{sub,1} = [actual_location,filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat'];
    pamlist{sub,2} = [actual_location,filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat'];

    if PAM_mirror_enabled == 1
    
        % here we can assign right state from the left, because fibers were
        % flipped, indices match
        pamlist{sub+numPatient,1} = [actual_location,filesep, 'stimulations',filesep,...
            ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat'];
        pamlist{sub+numPatient,2} = [actual_location,filesep, 'stimulations',filesep,...
            ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat'];
    end
end
