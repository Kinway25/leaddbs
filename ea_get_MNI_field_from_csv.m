function ea_get_MNI_field_from_csv(options, Field_array_file, Activation_threshold_VTA, sideLabel, templateOutputBasePath, source_index)

% converts any spatial field (e.g. VAT) to MNI based on available transformations
% the output resolution is defined by options.primarytemplate
% stored in [templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii']
% assumes OSS-DBS input and output

% VTA_array - 2D array, where the forth column is the field value
% Activation_threshold_VTA - V/m
% sideLabel - R or L 


% split to corrdinates and field
Field_array = table2array(readtable(Field_array_file, ReadVariableNames=false));
Field_coords = Field_array(:,2:4);
Field_vals = Field_array(:,8) * 1000.0;  % others are the components

% convert to native voxel space (will be as floating numbers)
Field_vox_native = ea_mm2vox(Field_coords, options.subj.preopAnat.(options.subj.AnchorModality).coreg)';

ea_dispt('Transforming spatial coordinates to MNI space...');

% apply normalization to MNI, not sure why we have to use the inverse
%Field_coords_MNI = ea_map_coords(Field_vox_native, options.subj.preopAnat.(options.subj.AnchorModality).coreg, ...
%            [options.subj.norm.transform.inverseBaseName,'ants.h5'], '', 'ANTS')';
Field_coords_MNI = ea_map_coords(Field_vox_native, options.subj.preopAnat.(options.subj.AnchorModality).coreg, ...
            [options.subj.subjDir, filesep, 'inverseTransform'], '')';



% now we need to deploy an equdistant grid in MNI space and interpolate across it

% just interpolate the magnitude, the vector field can be confusing
ea_dispt('Converting to equispaced image data...');
F = scatteredInterpolant(Field_coords_MNI(:,1),Field_coords_MNI(:,2),Field_coords_MNI(:,3),Field_vals,'linear','none');
gv=cell(3,1); spacing=zeros(3,1);

% hardwired N of points, if changed, also change Lattice shape in lead_settings.py
n_points = 71;   
for axis = 1:3
    %n_points(axis) = (max(round(Field_coords_MNI(:,axis))) - min(round(Field_coords_MNI(:,axis)))) / template.voxsize(axis);
    gv{axis}=linspace(min(round(Field_coords_MNI(:,axis))),max(round(Field_coords_MNI(:,axis))),n_points);
    spacing(axis)=abs(gv{axis}(1)-gv{axis}(2)); 
end

% I have no idea what is happening here
chun1=randperm(n_points); chun2=randperm(n_points); chun3=randperm(n_points); 
Vvat.mat=mldivide([(chun1);(chun2);(chun3);ones(1,n_points(1))]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,n_points)]')';

% My approach for MNI. Additional shift by half a voxel
% Vvat.mat = [template.voxsize(1), 0, 0, min(round(Field_coords_MNI(:,1)) - template.voxsize(1) / 2)
%             0, template.voxsize(2), 0, min(round(Field_coords_MNI(:,2)) - template.voxsize(2) / 2)
%             0, 0, template.voxsize(3), min(round(Field_coords_MNI(:,3)) - template.voxsize(3) / 2)
%             0, 0, 0, 1];

ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end
Vvat.dim=[n_points,n_points,n_points];
Vvat.dt = [4, endian];
Vvat.n=[1 1];
Vvat.descrip='oss-dbs-v2 - Field_ref';

ea_dispt('Filling data with values from interpolant...');
E_field_interp = F(gv);
E_field_interp(isnan(E_field_interp)) = 0;
E_field_interp(E_field_interp>10000.0) = 10000.0; % upperlimit files to 10000.
if source_index == 5  % no source indexing
    Vvat.fname = [templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    Vvat.fname = [templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_index), '.nii'];
end
Vvat.img = E_field_interp; 
ea_write_nii(Vvat);

% also create VATs directly
VTA_interp = E_field_interp;
VTA_interp = E_field_interp >= (Activation_threshold_VTA);
Vvat2 = Vvat;
Vvat2.descrip='oss-dbs-v2 - VAT_ref';
if source_index == 5  % no source indexing
    Vvat2.fname = [templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'];
else
    Vvat2.fname = [templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel,'_S',num2str(source_index), '.nii'];
end
Vvat2.pinfo = [1;0;352];
Vvat2.dt = [2, endian];
Vvat2.img = VTA_interp; 
ea_write_nii(Vvat2);
%ea_autocrop([templateOutputBasePath, 'binary_model-ossdbs_hemi-', sideLabel, '.nii'], margin=10);
%ea_autocrop([templateOutputBasePath, 'efield_model-ossdbs_hemi-', sideLabel, '.nii'], margin=10);




