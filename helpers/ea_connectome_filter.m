function ea_connectome_filter(connectome_path, ROI_file, threshold)
% This function allows to get fibers within the ROI
% By K.Butenko

arguments
    connectome_path                 % path to the folder containing the connectome files
    ROI_file % path to a nifti file containing a binary ROI that has to be intersected by the fibers to be preserved. Use a right hemisphere ROI, when flipping is used
    threshold % threshold to binarize ROI
end

% check if the downsampling parameter was provided
% and create an output folder

% sotre the filtered connectome next to the ROI
[ROI_folder,ROI_name,~] = fileparts(ROI_file);
C = strsplit(connectome_path,filesep);
connectome_name = C{end};
ROI_connectome_folder = [ROI_folder,filesep,connectome_name,'_',ROI_name];
mkdir(ROI_connectome_folder)

% load ROI to filter fibers
ROI = ea_load_nii(ROI_file);

%gets all mat files in struct
myFiles = dir(fullfile(connectome_path,'*.mat'));
% remove adjacency matrix and dataset_info if present
myFiles = myFiles(~endsWith({myFiles.name}, '_ADJ.mat'));
myFiles = myFiles(~endsWith({myFiles.name}, '_info.mat'));

for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    pathway_file = fullfile(myFiles(k).folder, baseFileName);

    new_pathway = [ROI_connectome_folder,filesep,myFiles(k).name];
    %fprintf(1, 'Now reading %s\n', fullFileName);

    ftr_full = load(pathway_file);

    % Trim connectome fibers by ROI
    ROI_Ind = find(abs(ROI.img(:))>threshold);   % ROI is assumed to be binary

    % Trim connectome fibers
    [xvox, yvox, zvox] = ind2sub(size(ROI.img), ROI_Ind);
    ROImm = ea_vox2mm([xvox, yvox, zvox], ROI.mat);
    filter = all(ftr_full.fibers(:,1:3)>=min(ROImm),2) & all(ftr_full.fibers(:,1:3)<=max(ROImm), 2);

    % discard the pathway if completely unconnected
    if ~any(filter)
        continue;
    end

    trimmedFiber = ftr_full.fibers(filter,:);

    [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
    fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, ROI)}, trimmedFiber(:,1:3), trimmedFiberID);
    
    % Remove outliers
    fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    connected = cellfun(@(fib) any(ismember(fib, ROI_Ind)), fibVoxInd);
    
    trimmedIdx = ftr_full.idx(trimmedFiberInd(connected),:);
    % restore complete trimmed fibers
    trimmedFiber = ftr_full.fibers(ismember(ftr_full.fibers(:,4), trimmedFiberInd(connected)), :);

    ftr.fibers = zeros(sum(trimmedIdx),4);
    ftr.idx = trimmedIdx;

    orig_indices = unique(trimmedFiber(:,4));
    for inx = 1:length(orig_indices)
        inx_to_change = trimmedFiber(:,4) == orig_indices(inx);
        trimmedFiber(inx_to_change,4) = inx;
    end

    ftr.fibers = trimmedFiber;

    % save the result
    ftr.fourindex = 1;
    ftr.ea_fibformat = '1.0';
    save(new_pathway, '-struct', 'ftr');

end

% % add mirror flag to already processed connectomes
% %gets all mat files in struct
% myFiles = dir(fullfile(connectome_path,'*.mat'));
% % remove adjacency matrix if present
% myFiles = myFiles(~endsWith({myFiles.name}, '_ADJ.mat'));
% 
% for k = 1:length(myFiles)
%     baseFileName = myFiles(k).name;
%     pathway_file = fullfile(myFiles(k).folder, baseFileName);
%     ftr = load(pathway_file);
% 
%     ftr.mirrored = 1;  % this flag designates that the connectome is fiberwise mirrorred
%     save(pathway_file, '-struct', 'ftr');
% end
