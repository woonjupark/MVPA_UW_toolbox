function [b] = VMPinVOI(vmp, voi, voiNum)
% [b] = GLMinVOI(vmp, voi ,voiNum)
%
% Inputs:
%   vmp             As given by vmp = BVQXfile(glmPath)
%   voi             As given by voi = BVQXfile(voiPath)
%   voiNum          Index number (or range) of the voi(s) to be indexed 
%                   (default: 1:length(voi.VOI))
% 
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