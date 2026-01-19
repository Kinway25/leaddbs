function [fib_index_in_WMH,fib_index_in_PVS,fib_index_in_lacunes] = ea_connectome_SVD_filter(ftr_full, WMH_file, PVS_file, lacunes_file)
% Find indices of fibers intersecting with different SVD
% By K.Butenko

arguments
    ftr_full                % (complete) ftr connectome (loaded)
    WMH_file % path to a nifti file containing a binary mask of WMH
    PVS_file % path to a nifti file containing a binary mask of WMH
    lacunes_file % path to a nifti file containing a binary mask of WMH
end


% process the fibers
if isfile(WMH_file)
    fib_index_in_WMH = trim_by_ROI(WMH_file,ftr_full);
else
    fib_index_in_WMH = false;
end

if isfile(PVS_file)
    fib_index_in_PVS = trim_by_ROI(PVS_file,ftr_full);
else
    fib_index_in_PVS = false;
end

if isfile(lacunes_file)
    fib_index_in_lacunes = trim_by_ROI(lacunes_file,ftr_full);
else
    fib_index_in_lacunes = false;
end

end

function fib_index_in_ROI = trim_by_ROI(ROI_file,ftr_full)
    threshold = 0.5; 
    ROI = ea_load_nii(ROI_file);
    
    % Trim connectome fibers by ROI
    ROI_Ind = find(abs(ROI.img(:))>threshold);   % ROI is assumed to be binary

    if isempty(ROI_Ind)
        fib_index_in_ROI = false;
        return
    end

    % Trim connectome fibers
    [xvox, yvox, zvox] = ind2sub(size(ROI.img), ROI_Ind);
    ROImm = ea_vox2mm([xvox, yvox, zvox], ROI.mat);
    filter = all(ftr_full.fibers(:,1:3)>=min(ROImm),2) & all(ftr_full.fibers(:,1:3)<=max(ROImm), 2);

    % discard the pathway if completely unconnected
    if ~any(filter)
        fib_index_in_ROI = false;
        return
    end

    trimmedFiber = ftr_full.fibers(filter,:);

    [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
    fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, ROI)}, trimmedFiber(:,1:3), trimmedFiberID);
    
    % Remove outliers
    fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    connected = cellfun(@(fib) any(ismember(fib, ROI_Ind)), fibVoxInd);
    
    fib_index_in_ROI = trimmedFiberInd(connected);
end
