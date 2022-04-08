function niiFiles = ea_dcm_to_nii(method, dicom_dir)

% simple wrapper that calls appropriate scripts to convert all dicoms in dicom_dir. Output is stored in tmp_dir
% method is an integer, 1 - 3
% 1 - dcm2nii, 2 - dicm2nii (Matlab), 3 - SPM)

tmp_dir = fullfile(dicom_dir, 'tmp');

% catch the case where method is a string
if ischar(method)
    switch method
        case 'dcm2niix'
            method_int = 1;
        case 'dicm2nii'
            method_int = 2;
        case 'SPM'
            method_int =3;
    end
else
    method_int = method;
end


switch method_int
    case 1 % dcm2niix
        ea_dcm2niix(dicom_dir, tmp_dir);
    case 2 % dicm2nii
        ea_dicm2nii(dicom_dir, tmp_dir);
    case 3 % SPM
        ea_spm_dicom_import(dicom_dir, tmp_dir);
end

% first option: rename files with the help of a GUI
[~, niiFiles] = fileparts(ea_regexpdir(tmp_dir, '\.nii\.gz$', 0));
if ischar(niiFiles)
    niiFiles = {niiFiles};
end

% clumsily remove .nii from filename
for idx = 1:length(niiFiles)
    niiFiles{idx, 1} = niiFiles{idx, 1}(1:end-4);
end

end