% mvpa.m
%
%  holds all support functions for mvpa toolbox
% project.
%
% functions can be called from outside with 'mvpa.<function name>'


classdef mvpa
    methods(Static)
        function roi = load_roi(paths)
            voiFile = fullfile(paths.roi{1}, [paths.subject, '_roi.voi']);
            roi = xff(voiFile); % loads in the roi file
        end
        function data_xff  = load_data(paths)
            disp(['loading ', paths]);
            clear data_xff;
            data_xff = xff(paths); % loading data into MATLAB

        end
        function emat = load_exp(fname, paths)
            scanInfo = regexp(fname, '.*_ses-(?<ses>\d+)_.*_run-(?<run>\d+)_.*', 'names');
            % load experimental protocols
            expPattern = sprintf('*_ses-%s_*_run-%s_*.mat', scanInfo.ses, scanInfo.run);
            expFile = dir(fullfile(paths.exp{1}, expPattern));
            load(fullfile(paths.exp{1}, expFile(1).name), 'emat');
        end
        function data_roi = subset_data_rois(data_xff, roi_xff, dataformat)
            % subset data with the rois
            if strcmp(dataformat, 'vmp')
                data_roi = mvpa.VMPinVOI(data_xff, roi_xff); % get the beta weights
            elseif strcmp(dataformat, 'glm')
                data_roi = mvpa.GLMinVOI(data_xff, roi_xff);% get the beta weights by roi
            end
        end
        function roi = collate_roi_predictors(roi, data_roi, dataformat)
            for r = 1:length(roi)
                roiIndx = strcmp({data_roi.name}, roi(r).name); % subset by the ROI
                if strcmp(dataformat, 'vmp')
                    roi(r).predictors = [roi(r).predictors; data_roi(roiIndx).beta'];  % block
                elseif strcmp(dataformat, 'glm')
                    roi(r).predictors = [roi(r).predictors; data_roi(roiIndx).beta(:,1:end-1)']; % block x timept
                end
            end
        end

        function factor = save_class_label(factor, emat)

            for f = 1:(length(factor)-1)
                factor(f).classlabels = cat(1, factor(f).classlabels, factor(f).labels(emat(:,factor(f).col))');
            end
                   
            for c = 1:size(emat,1)
                tmp = factor(1).classlabels{c};
                for f = 1:(length(factor)-1)
                    str{c} = strcat(tmp, factor(f).classlabels{c});
                end
            end
            factor(end).classlabels = cat(2, factor(end).classlabels, str);
        end

        function [b] = VMPinVOI(vmp, voi, voiNum)
            % [b] = GLMinVOI(vmp, voi ,voiNum)
            %
            % Inputs:
            %   vmp             As given by vmp = BVQXfile(glmpaths)
            %   voi             As given by voi = BVQXfile(voipaths)
            %   voiNum          Index number (or range) of the voi(s) to be indexed
            %                   (default: 1:length(voi.VOI))
            % Output:
            %   b              A structure containing vmp indices and data with fields:
            %       id         Linear indices of voi used for the vmp matrix
            %       beta       A [nVoxel nPredictors] matrix of beta weights within the
            %                  roi

            % Written by Kelly Chang - November 5, 2021

            %% Input Control

            if ~exist('voiNum', 'var')
                voiNum = 1:length(voi.VOI);
            end

            %% Start Extracting Beta Weights

            % vmp size and relative position (offset) within vmr volume
            map = vmp.Map; % extract map structure with vmp structure
            vmpData = cat(4, map.VMPData); % concatenate all vmp data

            % reshape vtc data into linear space
            vmpData = reshape(vmpData, [], size(vmpData, 4));

            vmpSize = size(vmp.Map(1).VMPData);
            vmpOffset = [vmp.XStart vmp.YStart vmp.ZStart];

            for i = voiNum
                % voi coordinates in anatomical resolution
                v = voi.BVCoords(i);

                % convert voi coordinates to vtc coordinates and resolution
                v = round(bsxfun(@minus, v, vmpOffset)/vmp.Resolution) + 1;

                % take only voi voxels inside the vtc volume
                indx = (v(:,1) > 0 & v(:,1) <= vmpSize(1) & ...
                    v(:,2) > 0 & v(:,2) <= vmpSize(2) & ...
                    v(:,3) > 0 & v(:,3) <= vmpSize(3));
                v = v(indx,:);

                % transform voi [x y z] coordinates into linear index equivalents
                v = sub2ind(vmpSize(1:3), v(:,1), v(:,2), v(:,3));
                v = unique(v); % keep only unique indices.

                % name of the roi
                b(i).name = voi.VOI(i).Name;

                % save the linear indices
                b(i).id = v;

                % record beta weight names
                b(i).conds = {map.Name};

                % take only vtc data inside voi voxels in linear space
                b(i).beta = vmpData(v, :);
            end
        end
        function [roi] = VTCinVOI(vtc, voi, voiNum, normalize)
            % [roi] = VTCinVOI(vtc, voi, voiNum, normalize)
            %
            % Extracts .VTC time courses within the specified .VOI region
            %
            % Input:
            %   vtc             As given by vtc = BVQXfile(vtcpaths)
            %   voi             As given by voi = BVQXfile(voipaths)
            %   voiNum          Index number (or range) of the voi(s) to be indexed
            %                   (default: 1:length(voi.VOI))
            %   normalize       Normalize vtc data: (vtc - mean(vtc)) / sd(vtc) (true)
            %                   OR not (false), logical (default: true)
            %
            % Output:
            %   roi             A structure containing VTC indices and data
            %       id          Linear indices of voi used for the vtc matrix
            %       vtcData     VTC data within the VOI (vtcData(:,roi.id))
            %
            % Note:
            % - Use [x,y,z] = ind2sub(vtcSize(2:end), roi.id) to transform from linear
            %   indices into 3D indices
            % - Dependencies: <a href="matlab: web('http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools/232-getting-started.html')">BVQXTools/NeuroElf</a>

            % Stolen primarily from vtc_VOITimeCourse from BVQXtools_v08b
            % Modified from getVOIcoordinates from Paola 20 Sep 2010
            % Edited by Kelly Chang - April 19, 2016

            %% Input Control

            % if voiNum does not exist, default = 1:length(voi.VOI)
            if ~exist('voiNum', 'var')
                voiNum = 1:length(voi.VOI);
            end

            % if normalize does not exist, default = true
            if ~exist('normalize', 'var')
                normalize = true;
            end

            %% Extract .vtc Data Located within .voi Space

            % vtc size and relative position (offset) within vmr volume
            vtcData = vtc.VTCData;
            vtcSize = size(vtcData);
            vtcOffset = [vtc.XStart vtc.YStart vtc.ZStart];

            % normalize vtc
            if normalize
                vtcData = bsxfun(@minus, vtcData, mean(vtcData)); % subtract each voxel by its mean
                vtcData = bsxfun(@rdivide, vtcData, std(vtcData)); % normalize each voxel by its SD
            end

            for i = voiNum
                % voi coordinates in anatomical resolution
                v = voi.BVCoords(i);

                % convert voi coordinates to vtc coordinates and resolution
                v = round(bsxfun(@minus, v, vtcOffset)/vtc.Resolution) + 1;

                % take only voi voxels inside the vtc volume
                indx = (v(:,1) > 0 & v(:,1) <= vtcSize(2) & ...
                    v(:,2) > 0 & v(:,2) <= vtcSize(3) & ...
                    v(:,3) > 0 & v(:,3) <= vtcSize(4));
                v = v(indx,:);

                % transform voi [x y z] coordinates into linear index equivalents
                v = sub2ind(vtcSize(2:end), v(:,1), v(:,2), v(:,3));

                % name of the roi
                roi(i).name = voi.VOI(i).Name;

                % only keep the unique indices
                roi(i).id = unique(v)';

                % reshape vtc data into linear space
                vtcData = reshape(vtcData, [vtcSize(1) prod(vtcSize(2:end))]);

                % take only vtc data inside voi voxels in linear space
                roi(i).vtcData = vtcData(:,roi(i).id);
            end
        end
        function [b] = GLMinVOI(glm, voi, voiNum)
            % [b] = GLMinVOI(glm, voi ,voiNum)
            %
            % Inputs:
            %   glm             As given by glm = BVQXfile(glmPath)
            %   voi             As given by voi = BVQXfile(voiPath)
            %   voiNum          Index number (or range) of the voi(s) to be indexed
            %                   (default: 1:length(voi.VOI))
            %
            % Output:
            %   b              A structure containing glm indices and data with fields:
            %       id         Linear indices of voi used for the glm matrix
            %       beta       A [nVoxel nPredictors] matrix of beta weights within the
            %                  roi

            % Written by Kelly Chang - October 10, 2016

            %% Input Control

            if ~exist('voiNum', 'var')
                voiNum = 1:length(voi.VOI);
            end

            %% Start Extracting Beta Weights

            % glm size and relative position (offset) within vmr volume
            glmData = glm.GLMData.BetaMaps;
            glmSize = size(glmData);
            glmOffset = [glm.XStart glm.YStart glm.ZStart];

            for i = voiNum
                % voi coordinates in anatomical resolution
                v = voi.BVCoords(i);

                % convert voi coordinates to vtc coordinates and resolution
                v = round(bsxfun(@minus, v, glmOffset)/glm.Resolution) + 1;

                % take only voi voxels inside the vtc volume
                indx = (v(:,1) > 0 & v(:,1) <= glmSize(1) & ...
                    v(:,2) > 0 & v(:,2) <= glmSize(2) & ...
                    v(:,3) > 0 & v(:,3) <= glmSize(3));
                v = v(indx,:);

                % transform voi [x y z] coordinates into linear index equivalents
                v = sub2ind(glmSize(1:3), v(:,1), v(:,2), v(:,3));

                % name of the roi
                b(i).name = voi.VOI(i).Name;

                % only keep the unique indices
                b(i).id = unique(v);

                % reshape vtc data into linear space
                glmData = reshape(glmData, [prod(glmSize(1:3)) glmSize(4)]);

                % take only vtc data inside voi voxels in linear space
                b(i).beta = glmData(b(i).id,:);
            end
        end
    end
end