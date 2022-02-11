clear all; close all;

% data: anatomical, functonal (task)
% roi: mt and pt
% experimental protocols (knowledge) which trial when
addpath(genpath('C:\Users\Ione Fine\Documents\code\neuroelf-matlab'))
%% directories

subject = 'sub-NS_G_RQ_1982';

paths.main = {'X:\WoonJuPark\MT\MTPilotTask'};
paths.roi = fullfile(paths.main, subject, 'ses-01', 'derivatives');
paths.beta = fullfile(paths.main, subject, {'ses-02', 'ses-03'}, 'derivatives');
paths.exp = fullfile(paths.main, subject);

%% load ROIs (aka .voi files)

voiFile = fullfile(paths.roi, [subject, '_roi.voi']);
roi = xff(voiFile); % HERE!!! this loads in the file

%% setup

% roi
rois = {'rPT'};

% conditions
condnames = {'Seq', 'Onset', 'Random'};
 
% conditions to include in training
includeConds = {'Seq', 'Onset', 'Random'};

% labels
labelid = 2; % 2: condition / 3: LR
if labelid == 2
    labelMatch = {'Seq', 'Onset', 'Random'}; % 1 = left, 2 = right    
elseif labelid == 3
    labelMatch = {'left', 'right'}; 
end

% beta from vmp or glm 
betaformat = 2; % 1: vmp; 2: glm

%% load beta weights

for whichroi = 1:length(rois)
    
    roiname = rois(whichroi); 

    designmat = []; % initialize the design matrix for roi voxels
    condsmat = {};  % initialize the conditions matrix for roi voxels
    labelsmat = {}; % initialize the labels matrix for roi voxels
    
    for i = 1:length(paths.beta) % for each session

        vmpFiles = dir(fullfile(paths.beta{i}, '*GLM_Trials.vmp'));
        glmFiles = dir(fullfile(paths.beta{i}, '*_glm_trials.glm'));

        for i2 = 1:length(vmpFiles) % for each vmp file

            % load vmp file
            if betaformat == 1
                
                vmpPath = fullfile(paths.beta{i}, vmpFiles(i2).name);
                vmp = xff(vmpPath); % loading vmp into MATLAB
                betaWeights = VMPinVOI(vmp, roi); % get the beta weights by roi
                
            elseif betaformat == 2
                
                glmPath = fullfile(paths.beta{i}, glmFiles(i2).name);
                glm = xff(glmPath); 
                betaWeights = GLMinVOI(glm, roi); 
                
            end

            % get session and run numbers to load exp protocols
            scanInfo = regexp(vmpFiles(i2).name, '.*_ses-(?<ses>\d+)_.*_run-(?<run>\d+)_.*', 'names');

            % load experimental protocols 
            expPattern = sprintf('*_ses-%s_*_run-%s_*.mat', scanInfo.ses, scanInfo.run);
            expFile = dir(fullfile(paths.exp{1}, expPattern));
            load(fullfile(paths.exp{1}, expFile(1).name), 'emat'); 

            % roi
            roiIndx = strcmp({betaWeights.name}, roiname);
            if betaformat == 1
                designmat = [designmat; betaWeights(roiIndx).beta'];
            elseif betaformat == 2
                designmat = [designmat; betaWeights(roiIndx).beta(:,1:end-1)'];
            end
            condsmat = cat(1, condsmat, condnames(emat(:,2))');
            labelsmat = cat(1, labelsmat, labelMatch(emat(:,labelid))');

        end
    end

    %% define train and test sets

    nconds = length(condnames);
    ntrialsPerRun = size(emat,1); 

    p = 1;

    for nruns = 2:12 % this is the cross-validation runs

        accuracy = nan(nruns, 3);
        ntrials = nan(nruns, 3);
        
        resultsMat = {{}, {}, {}};

        for testRun = 1:nruns % the fold

            % test set 
            testIndx = 1+(testRun-1)*ntrialsPerRun : ntrialsPerRun+(testRun-1)*ntrialsPerRun;

            testDesign = designmat(testIndx,:);
            testLabels = labelsmat(testIndx);
            testConds = condsmat(testIndx);

            % train set 
            allIndx = 1:ntrialsPerRun*(nruns-1)+ntrialsPerRun; 
            trainIndx = setdiff(allIndx, testIndx);

            tempConds = condsmat(trainIndx);
            tempDesign = designmat(trainIndx,:); 
            tempLabels = labelsmat(trainIndx); 

            condIndx = startsWith(tempConds, includeConds);

            trainDesign = tempDesign(condIndx,:);
            trainLabels = tempLabels(condIndx);
            trainConds = tempConds(condIndx);

            % run 
            [model,fitInfo] = fitcecoc(trainDesign, trainLabels); 
%             [model,fitInfo] = fitclinear(trainDesign, trainLabels); % linear classification model
            predictedLabels = predict(model, testDesign);

            results = [testLabels, predictedLabels];
                        
            % you can consider sem by trial sample (each run = 26 n's
            % add/removed)
            for i3 = 1:nconds

                thisIndx = startsWith(testConds, condnames(i3));
                accuracy(testRun,i3) = sum(strcmp(testLabels(thisIndx), predictedLabels(thisIndx)))/sum(thisIndx);
                ntrials(testRun,i3) = sum(thisIndx);
                
                resultsMat{i3} = cat(1, resultsMat{i3}, results(thisIndx,:));

            end

        end
        
        tmp = cellfun(@(x) strcmp(x(:,1), x(:,2)), resultsMat, 'UniformOutput', false);
        perfAvg(p,:) = cellfun(@mean, tmp); % averaging the accuracy value at the end
        perfSEM(p,:) = cellfun(@(x) std(x) / sqrt(length(x)), tmp); 

        p = p + 1;
        
    end

    figure(1);
    subplot(1,3,whichroi);
    errorbar(repmat([2:12]',1,3), perfAvg, perfSEM,'o-'); hold on; 
    plot(2:12, ones(1,11)*(1./length(labelMatch)), 'k--');
    legend(condnames);
    xlabel('nruns'); ylabel('accuracy');
    title(roiname);
    ylim([0 1]);
    axis('square');
    
end
