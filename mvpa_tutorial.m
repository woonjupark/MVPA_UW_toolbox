clear all; close all;

% Dependencies:
% Neuroelf toolbox
% https://github.com/neuroelf/neuroelf-matlab

addpath(genpath('C:\Users\Ione Fine\Documents\code\neuroelf-matlab'))

%% directories
paths.main = {'X:\WoonJuPark\MT\MTPilotTask'};
paths.subject = {'sub-NS_G_RQ_1982'};

%% load ROI (aka .voi files)
roi(1).name = {'rPT'}; roi(2).name = {'lMT'}; roi(3).name = {'rMT'};
paths.roi = {'rois'}; % where are the roi files located inside the subject directory
for run = 1:length(roi); roi(run).predictors = []; end % initialize the predictor matrix for roi voxels

%% load beta weights or glm data
paths.session = fullfile(paths.main, paths.subject, {'ses-02', 'ses-03'});

%% define if using vmp or glm BOLD data
% dataformat = 'vmp';
paths.data = fullfile('derivatives', '*_GLM_trials.vmp');
dataformat = 'glm';
paths.data = fullfile('derivatives', '*_glm_trials.glm'); % where the data files are located inside the subject directory

%% setup experimental condition lists
factor(1).col = 2;  factor(1).labels = {'Seq', 'Onset', 'Random'}; factor(1).chance = 1/3;
factor(2).col = 3; factor(2).labels =  {'left', 'right'}; factor(2).chance = 1/2;
factor(3).col = NaN; factor(3).labels = {'combo'}; factor(3).chance = 1/6; % combines the other factors
factor(4).col = NaN; factor(3).labels = {'session'}; factor(3).chance = NaN; % records session
factor(5).col = NaN; factor(4).labels = {'run'}; factor(3).chance = NaN; % records run

for f = 1:length(factor); factor(f).classlabels = [];  end

%% collect all the data
roi_xff = mvpa.load_roi(fullfile(paths.main, paths.subject,paths.roi,{[paths.subject{1}, '_roi.voi']})); % load the rois

for sess = 1:length(paths.session) % for each session
    cond_filelist = dir(fullfile(paths.session{sess}, '*task*.mat')); % each experimental condition file
    data_filelist = dir(fullfile(paths.session{sess}, paths.data)); % each data file
    if length(cond_filelist)~=length(data_filelist)
        error(['number of condition files ', num2str(length(condfilelist)), ...
            ' does not match the number of experimental files', num2str(length(data_filelist))]);
    end
    for run = 1:length(data_filelist) % for each vmp/glm file
        conds = mvpa.load_exp(cond_filelist(run)); % load exp protocols
        factor = mvpa.save_class_label(factor, conds, sess, run);        % save the class labels
        data_xff = mvpa.load_data(data_filelist(run)); % load in vmp or glm
        data_roi = mvpa.subset_data_rois(data_xff, roi_xff, dataformat); % just save the data for ROIs, in a temp structure
        roi = mvpa.collate_roi_predictors(roi, data_roi, dataformat); % now collate the roi data, over all the runs
    end
end

%% training time

model(1).desc = {'DiscrimType', 'linear', 'OptimizeHyperparameters', 'none'};
model(1).class_factor = 2; % which factor are you trying to classify
model(1).add_pred ={}; % {1, 'session', 'run'}; % factors, and other info to include as additional predictors
%       example models:
% 'DiscrimType: 'linear', 'quadratic', cubic
% 'OptimizeHyperparameters:,'none', 'auto};
%       if you want to do SVM this would be the terminology
% model(2).desc = {'SVM, OptimizeHyperparameters, none'};

model(1).CVstyle = {'Kfold', 10};
% other examples cross validation styles:
% {{'Kfold',  5}, {'Holdout', .1}, {'Leaveout', 'on'},
% if you want to generalize over specific conditions
% model(1).CVstyle= {'Generalize', 1, 'Seq', 'Onset'}; % define which factor are you Generalizing over, then the labels for train and sets

color_list = {'r', 'b', 'g' };
for r = 1:length(roi)
    for m = 1:length(model)
        predictors = mvpa.generate_predictors(model(m), factor, roi(r));
        if strcmp(model(m).CVstyle, 'Generalize')
            [perf, Mdl, Mdl_CV] = mvpa.classify(model(m),  roi(r).predictors, ...
                factor(model(m).class_factor).classlabels, factor(model(1).CVstyle{2}).classlabels);
        else
            [perf, Mdl, Mdl_CV] = mvpa.classify(model(m),  roi(r).predictors, ...
                factor(model(m).class_factor).classlabels);
        end
        h(m) = errorbar(r, perf.mean, perf.std, 'o'); 
        set(h(m), 'MarkerEdgeColor', color_list{m},'MarkerFaceColor', color_list{m}, 'Color', color_list{m}); hold on
    end
end
set(gca, 'XLim', [0 4])
