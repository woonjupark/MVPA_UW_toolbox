clear all; close all;

% Dependencies:
% Neuroelf toolbox
% https://github.com/neuroelf/neuroelf-matlab
addpath(genpath('C:\Users\Ione Fine\Documents\code\neuroelf-matlab'))
% data: anatomical, functonal (task)
% roi: mt and pt
% experimental protocols (knowledge) which trial when

%% directories
paths.main = {'X:\WoonJuPark\MT\MTPilotTask'};
paths.subject = 'sub-NS_G_RQ_1982';
%pathss.main = fullfile('~', 'Dropbox', '[WP]', '[Projects]', 'EB-MT', 'Data-MTPilotTask');

%% load ROI (aka .voi files)
roi(1).name = {'rPT'}; roi(2).name = {'lMT'}; roi(3).name = {'rMT'}; 
for i = 1:length(roi); roi(i).predictors = []; end % initialize the predictor matrix for roi voxels
paths.roi = fullfile(paths.main, paths.subject, 'ses-01', 'derivatives');
roi_xff = mvpa.load_roi(paths);

%% load beta weights or glm data
paths.session = fullfile(paths.main, paths.subject, {'ses-02', 'ses-03'}, 'derivatives');
% beta from vmp or glm
% dataformat = 'vmp';
% files.str = {'*GLM_Trials.vmp'};
dataformat = 'glm';
 files.str = {'*_glm_trials.glm'};

%% setup experimental condition lists
paths.exp = fullfile(paths.main, paths.subject);
% based on the experimental file's columns
factor(1).col = 2;  factor(1).labels = {'Seq', 'Onset', 'Random'}; factor(1).chance = 1/3;
factor(2).col = 3; factor(2).labels =  {'left', 'right'}; factor(2).chance = 1/2;
factor(3).col = NaN; factor(3).labels = NaN; factor(3).chance = 1/6;
for f = 1:length(factor); factor(f).classlabels = [];  end

% collect all the data
for s = 1:length(paths.session) % for each session
    datafilelist = dir(fullfile(paths.session{s}, files.str{1})); % each data file
    for i = 1:length(datafilelist) % for each vmp/glm file
        emat = mvpa.load_exp(datafilelist(i).name, paths); % load exp protocols
        % save the class labels
        factor = mvpa.save_class_label(factor, emat);
        data_xff = mvpa.load_data(fullfile(paths.session{s}, datafilelist(i).name)); % load in vmp or glm
        data_roi = mvpa.subset_data_rois(data_xff, roi_xff, dataformat); % just save the data for ROIs, in a temp structure
        roi = mvpa.collate_roi_predictors(roi, data_roi, dataformat); % now collate the roi data, over all the runs
    end
end

%% training time
model_list = {'fitcdiscr linear', 'fitcecoc'};
train_type ='generalize'; % 'kfold'; % generalize
if strcmp(train_type, 'kfold')
    train_idx = 1:length(factor(1).classlabels);
    test_idx = train_idx;
elseif strcmp(train_type, 'generalize')
    train_idx = strcmp(factor(1).classlabels, 'Seq');
    test_idx = strcmp(factor(1).classlabels, 'Onset');
end
color_list = {'r', 'b'};

for r = 1:length(roi)
    subplot(1, length(roi), r)
    for f = 1:length(factor)
        for m = 1:length(model_list)
            switch model_list{m}
                case 'fitcdiscr linear'
                    Mdl = fitcdiscr(roi(r).predictors(train_idx), factor(f).classlabels(train_idx),'OptimizeHyperparameters','none', 'DiscrimType', 'linear');
                case 'fitcecoc'
                    Mdl = fitcecoc(roi(r).predictors(train_idx), factor(f).classlabels(train_idx),'OptimizeHyperparameters','none');
            end
            %yfit = Mdl.predict(roi(r).predictors(test_idx));
            if strcmp(train_type, 'kfold')
                CVMdl = crossval(Mdl,  'CrossVal', 'on', 'KFold',10);
                perf(m,r,f).prob = 1-kfoldLoss(CVMdl,'Mode', 'individual');
                perf(m,r,f).mean = mean(perf(m,r,f).prob);
                perf(m,r,f).std = std(perf(m,r,f).prob ) / sqrt(length(perf(m,r,f).prob));
            else
                yfit = Mdl.predict(roi(r).predictors(test_idx));
                perf(m,r,f).mean = mean(strcmp(yfit, factor(2).classlabels(test_idx)));
                perf(m,r,f).std = 0 ;
            end
            h = errorbar(f, [perf(m,r,f).mean], [perf(m,r,f).std], 'o');
            set(h, 'MarkerEdgeColor', color_list{m},'MarkerFaceColor', color_list{m}); hold on
        end
        plot(f, factor(f).chance, 'k*');
    end
    set(gca, 'XLim', [0 length(factor)+1]);set(gca, 'XTick', 1:length(factor));
end



