function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, thresh)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

%pts_paths = ea_regexpdir('/media/konstantin/Konstantin/StimFit_Cohort/StimFitBIDS/derivatives/leaddbs','sub-*',0,'d',0);
%connectome_name = 'PetUpdatedMerged';

for side = 1:numSide

    pt_counter = 1;

%     if side == 1
%         side_tag = '_rh';
%     else
%         side_tag = '_lh';
%     end

    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalSum{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient

        if contains(vatlist{pt,side},'hemi-R')
            side_tag = '_rh';
        else
            side_tag = '_lh';
        end

        disp(['VAT ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if isstruct(vatlist) % direct nifti structs supplied
            vat = vatlist(pt,side);
        elseif iscell(vatlist) % filenames
            if isfile(vatlist{pt,side})
                vat = ea_load_nii(vatlist{pt,side});
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: VTA doesn''t exist!\n');
                continue;
            end
        end
        % Threshold the vat efield
        vatInd = find(abs(vat.img(:))>thresh);

        % Trim connectome fibers
        [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
        vatmm = ea_vox2mm([xvox, yvox, zvox], vat.mat);
        filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);

        % Skip further calculation in case VAT is totally not connected
        if ~any(filter)
            continue;
        end

        trimmedFiber = fibers(filter,:);

        % Map mm connectome fibers into VAT voxel space
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
        fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID);

        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

        % Find connected fibers
        connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);

        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(trimmedFiberInd(connected), pt)=1;

        % Checck intersection between vat and the connected fibers
        vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);

        % Generate fibsval for the Spearman's correlation method

        % instead of 5Peak, use max(E-proj)
        [~,pts,~] = fileparts(vatlist{pt,side}(1:end-35));
        pts_path = ['/media/konstantin/Konstantin/StimFit_Cohort/StimFitBIDS/derivatives/leaddbs/',pts];
        if contains(vatlist{pt,side},'20230506011354')
            file_E_proj = [pts_path,filesep,'miscellaneous',filesep,connectome_name,filesep,'Eproj_20230506011354',side_tag,filesep,'E_peak.mat'];
        else
            file_E_proj = [pts_path,filesep,'miscellaneous',filesep,connectome_name,filesep,'Eproj_06stimfit',side_tag,filesep,'E_peak.mat'];
        end
        E_proj = load(file_E_proj);
        fibsval5Peak{side}(trimmedFiberInd(connected), pt) = E_proj.E_peak(trimmedFiberInd(connected))*1000.0;
        pt_counter = pt_counter + 1;

        if rem(pt_counter-1,35) 
            pt_counter = 1;
        end

        fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);
        fibsvalPeak{side}(trimmedFiberInd(connected), pt) = cellfun(@max, vals);
        %fibsval5Peak{side}(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), vals);
    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    fibsvalSum{side} = sparse(fibsvalSum{side}(fibIsConnected, :));
    fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end
