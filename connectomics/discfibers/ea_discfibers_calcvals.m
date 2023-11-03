function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, thresh, map_list)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);


list_patients = [
"055fb1a649f7"
"1cJTEbdxV6EKuefRmS6bQ6a"
"22b959597732"
"23d9b814e2d9"
"2kZBXs9Tre3o9NAiQWQCb8"
"4d1f5973eb20"
"4eVfR2fjrSw2tMQLf9BbQt"
"8aPF46swc33pKnC6Rs4RaS"
"8c9c1101b527"
"8p3ECwNZ4WYyumErxvDmNx"
"91VmEGp8Z8x1ebK6N8kE1h"
"9EncFANdbNeJPU9p277oFH"
"9gq3BdcisR2T2qGvtjhf6Z"
"SBN0BBDF"
"SBN0DS64"
"SBN1S74E"
"SBN264D6"
"SBN2EB3F"
"SBN33DEA"
"SBN39922"
"SBN4E2BE"
"SBN6F007"
"SBN7F1C4"
"SBN9974A"
"SBN9F7S6"
"SBNA396B"
"SBNA94FS"
"SBNB0SBB"
"SBNB304A"
"SBNB40B6"
"SBNC46CB"
"SBNCB20B"
"SBND0360"
"SBNE1A7F"
"SBNEA29F"
"SBNEBCA3"
"SBNF72SA"
"SBNFSSFA"
"SBNS2366"
"b1y7LXWMe4dFu8nibrV8rL"
"b39f3f6874d3"
"beaucMTKeenmkGCfPCx4AQ"
"c9LM7G63AsUvxvv1T4cj8z"
"cFfPeBrdMafcN4TFCUExGE"
"cbiT8aFPJ67uFsTte1fAYw"
"d7FtVruwAebaqYM8CJKWFe"
"d8LCBe2YCXSJsc8QaTwQqE"
"eLQ1iuxG8nQDBr9yreUiGq"
"f1317979a99e"
"kY9M6btomZ9jELQQ8WA1wq"
"mvwcN7XHeFRkCkMeNVuupT"
"oKcxih6DX58fDpzHxzN1y3"
"ogqgsgAAsgKi2NAv3iqegW"
"sKFv3XzCKraSMEqnc4zRzJ"
"sMKbfUxmYtXmA2q6KxRveR"
"uWmYkD8MNLoRzAP7ydzh1R"
"uXSVxdRjFN5Mkmz3jbxyWV"
];


fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

for side = 1:numSide
    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalSum{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
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
        %fibsvalBin{side}(:, pt)=1;

        % Checck intersection between vat and the connected fibers
        vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);

        % Generate fibsval for the Spearman's correlation method
        fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);
        fibsvalPeak{side}(trimmedFiberInd(connected), pt) = cellfun(@max, vals);
        fibsval5Peak{side}(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), vals);

%         flag_mirror = 0;
% 
%         for pt_j = 1:size(list_patients,1)
%             % check which electrode
% 
%             temp2=strsplit(vatlist{pt},'_');
%             pt_ID = temp2{2}(5:end); 
% %             if contains(vatlist{i}, '_fl.')
% %                 AmpTrajDepthSide = [temp2{end-1}(4:end),'_lh'];
% %             else
% %                 AmpTrajDepthSide = [temp2{end}(4:end),'_lh'];
% %             end
% 
%             if strcmp(pt_ID, list_patients(pt_j)) 
% 
%                 AmptTrajectDepth = strsplit(vatlist{pt,side},'_');
%                 if contains(vatlist{pt,side}, 'right')
%                     AmptTrajectDepth = AmptTrajectDepth{end};
%                     file_E_proj = strcat(['/media/konstantin/Konstantin/JR/',char(list_patients(pt_j)),'/miscellaneous/PetUpdatedMerged/E_field_solution_',AmptTrajectDepth(1:end-4),'_rh/E_peak.mat']);
%                 else
%                     AmptTrajectDepth = AmptTrajectDepth{end-1};
%                     file_E_proj = strcat(['/media/konstantin/Konstantin/JR/',char(list_patients(pt_j)),'/miscellaneous/PetUpdatedMerged/E_field_solution_',AmptTrajectDepth,'_lh/E_peak.mat']);
%                     flag_mirror = 1;
%                 end
%                 E_proj = load(file_E_proj);
%                 break
%             end
%         end
%         %fibsval5Peak{side}(trimmedFiberInd(connected), pt) = E_proj.E_peak(trimmedFiberInd(connected));
% 
%             
%         if flag_mirror == 1
%             % for mirrored we load blocks of pathway counterparts as defined in
%             % obj.map_list (order is path1_rh,path1_lh,path2_rh...)
%             for pathway_i = 1:length(map_list)
%                 path_start = map_list(pathway_i);
% 
%                 if pathway_i ~= length(map_list)
%                     path_end = map_list(pathway_i+1) - 1;
%                 end
% 
%                 % odd numbers are rh, even lh
%                 if rem(pathway_i,2)
%                     path_start_counter = map_list(pathway_i+1);
%                     if pathway_i == length(map_list)-1
%                         disp("prelast pathway")
%                     else
%                         path_end_counter = map_list(pathway_i+2) - 1;
%                     end
%                 else
%                     path_start_counter = map_list(pathway_i-1);
%                     path_end_counter = map_list(pathway_i) - 1;                        
%                 end
% 
%                 % fiber statuses are stored in fibers, so we have check
%                 % first points on the fiber
%                 if pathway_i == length(map_list)-1
%                     fibsval5Peak{side}(path_start:path_end, pt) = E_proj.E_peak(path_start_counter:end);
%                 elseif pathway_i == length(map_list)
%                     fibsval5Peak{side}(path_start:end, pt) = E_proj.E_peak(path_start_counter:path_end_counter);
%                 else
%                     fibsval5Peak{side}(path_start:path_end, pt) = E_proj.E_peak(path_start_counter:path_end_counter);
%                 end
%                 %last_loc_i = fib_state_raw.idx(fib_i)+last_loc_i;            
%             end
%         else
%             fibsval5Peak{side}(:, pt) = E_proj.E_peak(:);
%         end

        %fibsvalBin{side}(:, pt) = (fibsval5Peak{side}(:, pt) >= 0.05);
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