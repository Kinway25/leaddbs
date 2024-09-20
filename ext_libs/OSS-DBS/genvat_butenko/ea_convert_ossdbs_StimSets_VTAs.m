function ea_convert_ossdbs_StimSets_VTAs(settings,side,outputPaths)
% Prepare Lead-DBS BIDS format VATs for unit contact-ground solutions.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    settings            % parameters for OSS-DBS simulation
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
end

switch side
    case 0
        sideLabel = 'R';
    case 1
        sideLabel = 'L';
end

% check how this works for only left electrodes
N_contacts = size(settings.contactLocation{1,side+1},1);

for contact_i = 1:N_contacts
    % also create 4D nii (4-th dimension is for E-field components and magnitude)
    file2save = [outputPaths.outputBasePath, '4D_efield_model-ossdbs_hemi-', sideLabel,'_',num2str(contact_i), '.nii'];
    ea_get_4Dfield_from_csv([outputPaths.HemiSimFolder, filesep, 'ResultsE1C',num2str(contact_i), filesep,'E_field_Lattice.csv'], file2save)
end