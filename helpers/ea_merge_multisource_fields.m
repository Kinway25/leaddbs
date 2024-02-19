function ea_merge_multisource_fields(basepath,source_efields,side,Activation_threshold_VTA,sideLabel)
  
%     copyfile([ea_space,'bb.nii'],[basepath,'bb_zero.nii']);
%     nii=ea_load_nii([basepath,'bb_zero.nii']);
%     nii.dt(1) = 16;
%     nii.img(:)=0;
%     ea_write_nii(nii);
%     
%     source_efields{side,end+1} = [basepath,'bb_zero.nii'];

    source_efields_side = {};
    cnt = 1;
    for i = 1:size(source_efields,2)
        % reverse order to add bb first
        if ~isempty(source_efields{side,end-i+1})
            source_efields_side{cnt,1} = source_efields{side,end-i+1};
            cnt = cnt + 1;
        end
    end

    output_dir = basepath;
    Vo = [basepath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'];
    matlabbatch{1}.spm.util.imcalc.input = source_efields_side;
    matlabbatch{1}.spm.util.imcalc.output = Vo;
    matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
    matlabbatch{1}.spm.util.imcalc.expression = 'max(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', matlabbatch);

    ea_delete([basepath,'bb_zero.nii'])

    % re-load, binarize and save
    nii=ea_load_nii(Vo);
    nii.img(:) = nii.img(:) >= Activation_threshold_VTA;
    nii.descrip='oss-dbs-v2 - VAT_ref';
    nii.fname = [basepath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
    ea_write_nii(nii);