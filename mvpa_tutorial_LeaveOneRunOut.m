% run MVPA using MVPA UW toolbox
% tutorial written by WJP 02/2022

clear all;

%% dependencies

% in addition to Neuroelf:
addpath('/Users/woonjupark/Dropbox/[WP]/[Projects]/Tools/MVPA_UW_Toolbox');

%% directories

paths.main = {'/Volumes/viscog.ibic.washington.edu/WoonJuPark/MT/MTPilotTask'};
paths.subject = {'sub-NS_G_RQ_1982'};

%% load ROI (aka .voi files)

roi(1).name = {'rPT'}; 
roi(2).name = {'lMT'}; 
roi(3).name = {'rMT'};

paths.roi = {'rois'}; % where are the roi files located inside the subject directory
for run = 1:length(roi)
    roi(run).predictors = [];
end % initialize the roi struct to collate later 

%% load beta weights

paths.session = fullfile(paths.main, paths.subject, {'ses-02', 'ses-03'});

%% define if using vmp or glm BOLD data

% dataformat = 'vmp';
% paths.data = fullfile('derivatives', '*_GLM_trials.vmp');
dataformat = 'glm';
paths.data = fullfile('derivatives', '*_glm_trials_ver2.glm'); % where the data files are located inside the subject directory

%% setup experimental condition lists
% make sure to add 'session' and 'run', as some of the info are used in other functions

factor(1).col = 2; 
factor(1).labels = {'Seq', 'Onset', 'Random'}; 
factor(1).chance = 1/3;

factor(2).col = 3; 
factor(2).labels =  {'left', 'right'}; 
factor(2).chance = 1/2;

factor(3).col = NaN; 
factor(3).labels = {'combo'}; 
factor(3).chance = 1/6; % combines the other factors

factor(4).col = NaN; 
factor(4).labels = {'session'}; 
factor(4).chance = NaN; % records session

factor(5).col = NaN; 
factor(5).labels = {'run'}; 
factor(5).chance = NaN; % records run

for f = 1:length(factor)
    factor(f).classlabels = [];
end % initialize the factors to collate later 

%% collect all the rois

voiFile = fullfile(paths.main, paths.subject,paths.roi,{[paths.subject{1}, '_undist_rois.voi']});
roi_xff = mvpa.load_roi(voiFile{1}); % load the rois

for sess = 1:length(paths.session) % for each session
    cond_filelist = dir(fullfile(paths.session{sess},'*task*.mat')); % each experimental condition file
    data_filelist = dir(fullfile(paths.session{sess}, paths.data)); % each data file
    if length(cond_filelist)~=length(data_filelist)
        error(['number of condition files ', num2str(length(cond_filelist)), ...
            ' does not match the number of experimental files', num2str(length(data_filelist))]);
    end
    for run = 1:length(data_filelist) % for each vmp/glm file

        % deal with factors 
        conds = mvpa.load_exp(cond_filelist(run)); % load exp protocols
        factor = mvpa.collate_factor_labels(factor, cell2mat(conds.mat), sess, run); % save the class labels and make everything lowercase
        
        % deal with data
        data_xff = mvpa.load_data(data_filelist(run)); % load in vmp or glm
        data_roi = mvpa.subset_data_rois(data_xff, roi_xff, dataformat); % just save the data for ROIs, in a temp structure
        roi = mvpa.collate_roi_predictors(roi, data_roi, dataformat); % now collate the roi data, over all the runs
    end
end


%% classify L/R using Leave One Run Out procedure for cross validation

model(1).desc = {'DiscrimType', 'linear', 'OptimizeHyperparameters', 'none'};
model(1).class_factor = 2; % which factor are you trying to classify? (1: stim type; 2: L/R)
model(1).add_pred = {'session','run'};
model(1).CVstyle= {'LeaveOneRunOut'}; 
model(1).trainConds = {'seq'};
model(1).testConds = {'seq','onset','random'};
model(1).cond_factor = 1; 
model(1).color = 'r'; 
model(1).sym = 's';
model(1).figname = 'Classify L/R, train: seq';

model(2).desc = {'DiscrimType', 'linear', 'OptimizeHyperparameters', 'none'};
model(2).class_factor = 1; % which factor are you trying to classify? (1: stim type; 2: L/R)
model(2).add_pred = {'session','run'};
model(2).CVstyle= {'LeaveOneRunOut'}; 
% model(2).CVstyle= {'Kfold',10}; % can do Kfold if train and test conds are the same
model(2).trainConds = {'seq','onset','random'};
model(2).testConds = {'seq','onset','random'};
model(2).cond_factor = 1; 
model(2).color = 'b'; 
model(2).sym = 'o';
x_lim = [0 4];
x_chance = min(x_lim):max(x_lim);
model(2).figname = 'Classify StimType, train: all';

for m = 1:2
    
    figure('NumberTitle', 'off', 'Name', model(m).figname);
    
    for r = 1:3
        
        % generate predictors 
        % if not adding factors, then this should work: predictors = roi(r).predictors;         
        predictors = mvpa.generate_predictors(model(m), factor, roi(r)); % add additional factors (e.g., run/sessions)
        predictors = cell2mat(predictors);
        
        % partition data for cross validation
        [cv, doi, doiPredictors] = mvpa.select_train_test(model(m), factor, predictors);
        
        % initialize
        testlabels_all = [];
        predlabels_all = [];
        coi_all = [];
        
        % train and test 
        for run = 1:cv.NumTestSets % run will be the test set
            
            % get train and test sets
            train_idx = cv.training(run);
            test_idx = cv.test(run);
            trainlabels = doi(model(m).class_factor).doilabels(train_idx);
            testlabels = doi(model(m).class_factor).doilabels(test_idx);
            coi = doi(model(m).cond_factor).doilabels(test_idx);
            
            % train
            Mdl = mvpa.train(model(m), doiPredictors(train_idx,:), trainlabels);
            
            % test
            predlabels = mvpa.test(Mdl, doiPredictors(test_idx,:));
            testlabels_all = [testlabels_all; testlabels];
            predlabels_all = [predlabels_all; predlabels];
            coi_all = [coi_all; coi];
        end
        
        % get performance and plot
        for c = 1:length(model(m).testConds)
            condidx = strcmp(coi_all,model(m).testConds{c});
            results(c) = mvpa.output(testlabels_all(condidx), predlabels_all(condidx));
            subplot(1,3,r);
            h = errorbar(c, results(c).mean, results(c).sem, model(m).sym); hold on;
            set(h, 'MarkerEdgeColor', model(m).color, 'MarkerFaceColor', model(m).color, 'Color', model(m).color); 
        end
        plot(x_chance, repmat(factor(model(m).class_factor).chance,size(x_chance)), 'k--');
        set(gca, 'XLim', x_lim)
        set(gca, 'YLim', [0 1])
        title(roi(r).name);
        xlabel('conditions')
    end
end



