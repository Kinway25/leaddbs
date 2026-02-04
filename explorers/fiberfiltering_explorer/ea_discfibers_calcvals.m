function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, thresh)
% Calculate fiber connection values based on the VATs and the connectome
% WARNING: modified for SVD flipped to one hemisphere (numSide=1)

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);

if numSide ~= 1
    ea_warndlg("The SVD script expects single stimulation per protocol!!!")
    return
end

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

WMH_lacunes_fibers = zeros(numPatient,2); % total and stim
PVS_fibers = zeros(numPatient,2);

for side = 1:numSide
    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalSum{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    % Because PVS are treated as NaNs, we need to include them for
    % non-connected cases (if this fiber is connected to another stims)
    % to avoid treating them as 0s!
    fibsvalPVS{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
        disp(['VAT ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if isstruct(vatlist) % direct nifti structs supplied
            vat = vatlist(pt,side);
        elseif iscell(vatlist) % filenames
            if strcmp(vatlist{pt,side},"skip")
                % no stimulation for this hemisphere
                continue
            elseif isfile(vatlist{pt,side})
                vat = ea_load_nii(vatlist{pt,side});
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: VTA doesn''t exist!\n');
                continue;
            end
        end
        % Threshold the vat efield
        vatInd = find(abs(vat.img(:))>thresh);

        % bb including SVDs 
        vatInd_ext = find(abs(vat.img(:))>thresh | vat.img(:) == -1 | vat.img(:) == -2);
        %vatInd_ext = find(abs(vat.img(:))>thresh);

        % Trim connectome fibers
        %[xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
        [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd_ext);
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

        % Checck intersection between vat and the connected fibers
        %vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);
        vals = cellfun(@(fib) vat.img(intersect(fib, vatInd_ext)), fibVoxInd(connected), 'Uni', 0);

        % count all fibers
        connected_SVD = cellfun(@(fib) any(ismember(fib, vatInd_ext)), fibVoxInd);
        vals_all = cellfun(@(fib) vat.img(intersect(fib, vatInd_ext)), fibVoxInd(connected_SVD), 'Uni', 0);
        vals_all_min = cellfun(@min, vals_all);
        % WMH_lacunes_fibers(pt,1) = sum(vals_all_min == -1);
        % PVS_fibers(pt,1) = sum(vals_all_min == -2);
        trimmedFiberInd_connSVD = trimmedFiberInd(connected_SVD);
        fibsvalPVS{side}(trimmedFiberInd_connSVD(vals_all_min == -2),pt) = 1;


        % SVD correction
        vals_min = cellfun(@min, vals);

        fibsvalBin{side}(trimmedFiberInd(connected), pt)=1;
        fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);
        fibsvalPeak{side}(trimmedFiberInd(connected), pt) = cellfun(@max, vals);
        fibsval5Peak{side}(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), vals);

        trimmedFiberInd_conn = trimmedFiberInd(connected);

        if any(vals_min == -1)

            %WMH_lacunes_fibers(pt,2) = sum(vals_min == -1);

            disp("WMH/lacunes intersection detected")
            fibsvalBin{side}(trimmedFiberInd_conn(vals_min < -0.9), pt) = 0;
            fibsvalSum{side}(trimmedFiberInd_conn(vals_min < -0.9), pt) = 0;
            fibsvalMean{side}(trimmedFiberInd_conn(vals_min < -0.9), pt) = 0;
            fibsvalPeak{side}(trimmedFiberInd_conn(vals_min < -0.9), pt) = 0;
            fibsval5Peak{side}(trimmedFiberInd_conn(vals_min < -0.9), pt) = 0;
        end

        if any(vals_min == -2)

            %PVS_fibers(pt,2) = sum(vals_min == -2);

            disp("PVS intersection detected")
            fibsvalBin{side}(trimmedFiberInd_conn(vals_min == -2), pt) = nan;
            fibsvalSum{side}(trimmedFiberInd_conn(vals_min == -2), pt) = nan;
            fibsvalMean{side}(trimmedFiberInd_conn(vals_min == -2), pt) = nan;
            fibsvalPeak{side}(trimmedFiberInd_conn(vals_min == -2), pt) = nan;
            fibsval5Peak{side}(trimmedFiberInd_conn(vals_min == -2), pt) = nan;
        end
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

    % add PVS NaNs for fibers that were connected in other stims
    fibsvalBin{side}(logical(fibsvalPVS{side}(fibIsConnected,:))) = nan;
    fibsvalSum{side}(logical(fibsvalPVS{side}(fibIsConnected,:))) = nan;
    fibsvalMean{side}(logical(fibsvalPVS{side}(fibIsConnected,:))) = nan;
    fibsvalPeak{side}(logical(fibsvalPVS{side}(fibIsConnected,:))) = nan;
    fibsval5Peak{side}(logical(fibsvalPVS{side}(fibIsConnected,:))) = nan;

end
