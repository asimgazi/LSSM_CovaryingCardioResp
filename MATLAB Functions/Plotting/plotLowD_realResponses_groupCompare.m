function plotLowD_realResponses_groupCompare(dim2plot, plotTitle, fileName, ...
    featSegm_groupSep, Fs, varargin)
% This function will take real responses of different groups of data,
% stack them all up, apply PCA, and then reseparate according to groups
% according to inputs provided

% ---------- Inputs -------------------- %
% dim2plot: scalar = 2 or 3 that dictates whether to plot in 2-D or 3-D
% plotTitle: figure title to use
% fileName: filename to use
% featSegm_groupSep: cell array of cell arrays. First layer of cells
%                    separates groups. Second layer of cells separates
%                    participants within the groups. Each cell contains a
%                    struct that contains a field for each protocol
%                    condition, where within each field, an array exists
%                    that is TxMxC, where T is the number of time
%                    points in that condition segment, M is the number of
%                    features, and C is the number of condition repetitions
%                    for that condition
% Fs: sampling frequency; used to construct time vectors to plot against

% (optional)
% 'plotOption': flag that specifies to look for the option for plotting in
%               next function input. Options to choose from:
%               'condSep'
%               'condAvg'
%               'condAvg_train': if this option is selected, then the next
%               input should be the train-test split condition code
%               (e.g., 6.5 for split between 6 and 7)

% 'subjOverlay' - if this flag is provided, the next varargin cell needs to
% be a cell array with a vector of subject numbers for each group in each
% cell to use for legends. What this flag entails is it will set
% meanResp = false. This will overlay all subjects' responses from a group for a
% particular condition on top of each other while separating groups and
% conditions into separate plots. If meanResp = true, means will be computed across
% participants for each group and condition and then overlaid on the same plot for
% easy comparison, so the legend will just be based on condition and group

% 'verbose' if you want to see all warnings (default is not verbose)

% 'groupNames' will need to be followed by a cell array of groups names
% corresponding to responses_groupSep


% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'condSep')
        condAvg = false;
        trainOnly = false;
    elseif strcmp(varargin{arg}, 'condAvg')
        condAvg = true;
        trainOnly = false;
    elseif strcmp(varargin{arg}, 'condAvg_train')
        condAvg = true;
        trainOnly = true;
        trTest_split = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'verbose')
        verbose = true;
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanResp = false;    % We will plot responses themselves
        subjectNums = varargin{arg + 1};    % Subject numbers for legend
    elseif strcmp(varargin{arg}, 'groupNames')
        groupNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'condNames')
        condNames_legend = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('condAvg', 'var')
    condAvg = true;
end
if ~exist('trainOnly', 'var')
    trainOnly = true;
end
if ~exist('trTest_split', 'var')
    trTest_split = 6.5;
end
if ~exist('verbose', 'var')
    verbose = false;
end
if ~exist('meanResp', 'var')
    meanResp = true; % read documentation above for what this bool dictates
end
if ~exist('condNames_legend', 'var')
    condNames_legend = fieldnames(featSegm_groupSep{1}{1});
end

% We suppress these expected warnings if verbosity is not desired
% (we already ignore these models when optimizing HPs, so we might as well
% ignore the warnings)
if ~verbose
    % Specify the warning IDs to turn off
    wrnID = {'stats:pca:ColRankDefX'};
    
    % Turn these warnings off
    for id = 1:numel(wrnID)
        warning('off', wrnID{id})
    end
end


% Just to have these saved as easy to understand variables
numGroups = numel(featSegm_groupSep);
numSubjects = zeros(numGroups, 1);  % # subjects can differ between groups
for grp_idx = 1:numGroups
    numSubjects(grp_idx) = numel(featSegm_groupSep{grp_idx});
end
condNames = fieldnames(featSegm_groupSep{1}{1});
numConds = numel(condNames);
numSamples = zeros(numConds, 1);    % # samples can differ between conditions
for cond_idx = 1:numConds
    numSamples(cond_idx) = ...
        size(featSegm_groupSep{1}{1}.(condNames{cond_idx}), 1);
end
numFeats = size(featSegm_groupSep{1}{1}.(condNames{1}), 2);


% Branch based on whether we average within a subject across conditions
if condAvg
    % Initialize a data structure to hold condition-averaged responses for
    % each subject
    feat_condAvg_groupSep = cell(numGroups, 1);
    for grp = 1:numGroups
        % Groups can have different numbers of subjects
        feat_condAvg_groupSep{grp} = cell(numSubjects(grp), 1);
        for sub = 1:numSubjects(grp)
            feat_condAvg_groupSep{grp}{sub} = struct();
            for cond = 1:numConds
                % Initialized to NaNs as default in case no data are available
                % to replace (will ignore NaN responses during plotting)
                feat_condAvg_groupSep{grp}{sub}.(condNames{cond}) = ...
                    nan*ones(size(featSegm_groupSep{grp}{sub}.(condNames{cond}), 1), ...
                    numFeats);
            end
        end
    end
    
    % Depending on if we want to only use training data or not, we need to
    % average differently
    if trainOnly
        % This function is created for the DARPA PTSD dataset, so the
        % trTest_split convention follows the following:
        % |N1|N2|T1V1|T2V2|V3|V4|N3|N4|T3V5|T4V6| (10 conditions in total)
        % So, if a user inputs trTest_split == 6.5, the split will occur
        % after V4 before N3, for example
        if trTest_split > 3 && trTest_split < 4
            % Number of training repetitions to include for each condition
            trainV = 1;
            trainT = 1;
            trainN = 2;
        elseif trTest_split > 4 && trTest_split < 5
            trainV = 2;
            trainT = 2;
            trainN = 2;
        elseif trTest_split > 5 && trTest_split < 6
            trainV = 3;
            trainT = 2;
            trainN = 2;
        elseif trTest_split > 6 && trTest_split < 7
            trainV = 4;
            trainT = 2;
            trainN = 2;
        elseif trTest_split > 7 && trTest_split < 8
            trainV = 4;
            trainT = 2;
            trainN = 3;
        elseif trTest_split > 8 && trTest_split < 9
            trainV = 4;
            trainT = 2;
            trainN = 4;
        elseif trTest_split > 9 && trTest_split < 10
            trainV = 5;
            trainT = 3;
            trainN = 4;
        else
            error('This train-test split is not supported.');
        end
        
        % Now that we know how many V, T, and N to include, we can average
        trainNums = [trainV, trainT, trainN];
        for grp = 1:numGroups
            for sub = 1:numSubjects(grp)
                for cond = 1:numConds
                    % Average along the 3rd dimension to average over repetitions
                    % Only average the training data
                    feat_condAvg_groupSep{grp}{sub}.(condNames{cond}) = ...
                        mean(featSegm_groupSep{grp}{sub}.(condNames{cond})(:, :, ...
                        1:trainNums(cond)), 3, 'omitnan');
                end
            end
        end
    else
        % In this case, we can just average over whatever we have since no
        % train/test considerations
        for grp = 1:numGroups
            for sub = 1:numSubjects(grp)
                for cond = 1:numConds
                    % Average along the 3rd dimension to average over repetitions
                    feat_condAvg_groupSep{grp}{sub}.(condNames{cond}) = ...
                        mean(featSegm_groupSep{grp}{sub}.(condNames{cond}), 3, 'omitnan');
                end
            end
        end
    end
    
    % Separate responses by group and by condition; within each separation,
    % stack across all subjects
    resps_groupSep_condSep = cell(numGroups, numConds);
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            % Initialize
            resps_groupSep_condSep{grp_idx, cond_idx} = ...
                zeros(numSamples(cond_idx), numFeats, numSubjects(grp_idx));
        end
    end
    
    % Now fill the cell array's arrays
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            for sub = 1:numSubjects(grp_idx)
                resps_groupSep_condSep{grp_idx, cond_idx}(:, :, sub) = ...
                    feat_condAvg_groupSep{grp_idx}{sub}.(condNames{cond_idx});
            end
        end
    end
    
    % To encode time, I am going to create a time vector according to the
    % longest time protocol condition we have and then just use a subset of
    % that for the other conditions
    timeVec = 0:(1/Fs):(max(numSamples)-1)/Fs;
    
    % Encode time in both alpha and size of scatter plot points using
    % alpha and size conversion rule below
    alphaMin = 0.3; sizeMin = 20; timeMin = min(timeVec);
    alphaMax = 1; sizeMax = 70; timeMax = max(timeVec);
    
    % Change in transparency and size per time increment
    alphaSlope = (alphaMax - alphaMin)/(timeMax - timeMin);
    sizeSlope = (sizeMax - sizeMin)/(timeMax - timeMin);
    
    % Store away alpha and size "labels" to encode time
    alphaLabels = alphaSlope*(timeVec - timeMin) + alphaMin;
    sizeLabels = sizeSlope*(timeVec - timeMin) + sizeMin;
    
    % Branch based on meanResp
    if meanResp
        % Create basic color labels (blue, red, green) for protocol conditions
        colorLabels = (1/255)*[221, 85, 255; ...
            255, 85, 85; ...
            95, 211, 188];
        
        % To differentiate groups, we will use different markers
        grpLabels = {'o', 'p', 's', '^'};
        
        % Compute means within each condition, for each group
        % (3rd dimension is subject)
        meanResps_groupSep_condAvg = cell(numGroups, numConds);
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                meanResps_groupSep_condAvg{grp_idx}{cond_idx} = ...
                    mean(resps_groupSep_condSep{grp_idx, cond_idx}, 3);
            end
        end
        
        % Initialize stacked responses
        stackedResps = [];
        groupLabels = [];   % To easily re-separate groups
        condLabels = [];    % To easily re-separate conditions
        
        % Loop and stack
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                stackedResps = [stackedResps; ...
                    meanResps_groupSep_condAvg{grp_idx}{cond_idx}];
                groupLabels = [groupLabels; grp_idx*...
                    ones(size(meanResps_groupSep_condAvg{grp_idx}{cond_idx}, 1), 1)];
                condLabels = [condLabels; cond_idx*...
                    ones(size(meanResps_groupSep_condAvg{grp_idx}{cond_idx}, 1), 1)];
            end
        end
        
        % Check to make sure sizes are correct
        if size(stackedResps, 1) ~= sum(numSamples)*numGroups || ...
                length(condLabels) ~= sum(numSamples)*numGroups || ...
                length(groupLabels) ~= sum(numSamples)*numGroups
            error('Stacking did not go as planned');
        end
        
        % Now I can make each column unit variance for PCA
        pcaReady = zeros(size(stackedResps));
        for col = 1:size(stackedResps, 2)
            pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
        end
        
        % Apply PCA and store away variance explained and PCs
        % (outputs I ignore are the principal component coefficients, the variances
        % in the latent space, Hotelling's T-squared statistic for each
        % observation, and the mean of each column)
        [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % Separate the PCA scores by group and condition to make labeling of
        % plots easy (i.e., "unstack")
        plotReady = cell(numGroups, numConds);
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                plotReady{grp_idx, cond_idx} = ...
                    pcaResult(groupLabels == grp_idx & ...
                    condLabels == cond_idx, :);
                
                % Check size to make sure everything was stored properly
                if size(plotReady{grp_idx, cond_idx}, 1) ~= numSamples(cond_idx) || ...
                        size(plotReady{grp_idx, cond_idx}, 2) ~= numFeats
                    error('plotReady sizing is incorrect')
                end
            end
        end
        
        % Branch based on dimensions to plot
        if dim2plot == 2
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, numGroups);
            
            % Iterate through groups and conditions and overlay
            for grp_idx = 1:numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:numConds
                    % Plot one scatter point at a time to leverage time labels
                    for t = 1:size(plotReady{grp_idx, cond_idx}, 1)
                        s = scatter(plotReady{grp_idx, cond_idx}(t, 1), ...
                            plotReady{grp_idx, cond_idx}(t, 2), sizeLabels(t), ...
                            colorLabels(cond_idx, :), grpLabels{grp_idx}, 'filled');
                        s.MarkerFaceAlpha = alphaLabels(t);
                        
                        % If this is the final datapoint, then store legend
                        % info
                        if t == size(plotReady{grp_idx, cond_idx}, 1)
                            s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                        else
                            s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                
                % Include legend
                legend(condNames_legend, 'Location', 'best')
                legend('boxoff')
                hold off
                
                % Add title to subplot
                title(groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Create video for all data overlaid
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames_legend, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
            
            % Close all figures
            close all
            
            
        else
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, numGroups);
            
            % Iterate through groups and conditions and overlay
            for grp_idx = 1:numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:numConds
                    % Plot one scatter point at a time to leverage time labels
                    for t = 1:size(plotReady{grp_idx, cond_idx}, 1)
                        s = scatter3(plotReady{grp_idx, cond_idx}(t, 1), ...
                            plotReady{grp_idx, cond_idx}(t, 2), ...
                            plotReady{grp_idx, cond_idx}(t, 3), ...
                            sizeLabels(t), ...
                            colorLabels(cond_idx, :), grpLabels{grp_idx}, 'filled');
                        s.MarkerFaceAlpha = alphaLabels(t);
                        
                        % If this is the final datapoint, then store legend
                        % info
                        if t == size(plotReady{grp_idx, cond_idx}, 1)
                            s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                        else
                            s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                
                % Include legend
                legend(condNames_legend, 'Location', 'best')
                legend('boxoff')
                hold off
                
                % Add title to subplot
                title(groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Create video for all data overlaid
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames_legend, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
            
            % Close all figures
            close all
        end
        
        % Create bar graphs for each of the PCs plotted specifying their
        % "loadings" a.k.a. "recipes"
        fig = figure(1);
        set(fig, 'Visible', 'on');
        for dim = 1:dim2plot
            % Create a dim2plot x 1 subplot grid
            subplot(dim2plot, 1, dim)
            
            % Each principal component's loadings is a column of the matrix
            bar(pcaLoads(:, dim), 'red')
            
            % Label axes and title
            xlabel('Feature Index')
            ylabel('Coefficient')
            title(['PC', num2str(dim), ' Recipe'])
        end
        % Set title for entire figure
        sgtitle('Principal Component Loadings')
        
        % Save this figure away
        saveas(fig, [fileName, '_PCloadings'], 'fig');
        
        % Close all plots
        close all
        
    else
        % Create color map based on total number of subjects
        rng(0);     % Set random seed so that colors are always the same
        colorLabels = rand(sum(numSubjects), 3);
        
        % To differentiate groups, we will use different markers
        grpLabels = {'o', 'p', 's', '^'};
        
        % Create nested cell array of subject IDs based on groups and subject
        % numbers
        subjectIDs = cell(numGroups, 1);
        for grp_idx = 1:numGroups
            subjectIDs{grp_idx} = cell(numSubjects(grp_idx), 1);
        end
        for grp_idx = 1:numGroups
            for sub = 1:numel(subjectIDs{grp_idx})
                subjectIDs{grp_idx}{sub} = num2str(subjectNums{grp_idx}(sub));
            end
        end
        
        % Since I am not doing any averaging across subjects of the same
        % group, I essentially need to stack and flatten all of the
        % responses into a 2-D array of dimension:
        % sum(numSubjects)*sum(numSamples) x numFeats
        % At the same time, I will need to keep track of groupLabels,
        % subjLabels, and condLabels so I can later separate when plotting
        groupLabels = [];
        subjLabels = [];
        condLabels = [];
        
        % Initialize the stacked result
        stackedResps = [];
        
        % Loop and stack
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                for sub_idx = 1:numSubjects(grp_idx)
                    % Stack responses
                    stackedResps = [stackedResps; ...
                        resps_groupSep_condSep{grp_idx, cond_idx}(:, :, sub_idx)];
                    
                    % Keep track of grouping labels
                    groupLabels = [groupLabels; grp_idx*...
                        ones(size(resps_groupSep_condSep{grp_idx, cond_idx}...
                        (:, :, sub_idx), 1), 1)];
                    subjLabels = [subjLabels; subjectNums{grp_idx}(sub_idx)*...
                        ones(size(resps_groupSep_condSep{grp_idx, cond_idx}...
                        (:, :, sub_idx), 1), 1)];
                    condLabels = [condLabels; cond_idx*...
                        ones(size(resps_groupSep_condSep{grp_idx, cond_idx}...
                        (:, :, sub_idx), 1), 1)];
                end
            end
        end
        
        % Check to make sure sizes are correct
        if size(stackedResps, 1) ~= sum(numSubjects)*sum(numSamples) || ...
                length(subjLabels) ~=  sum(numSubjects)*sum(numSamples) || ...
                length(condLabels) ~= sum(numSubjects)*sum(numSamples)
            error('Stacking did not go as planned');
        end
        
        % Now I can make each column unit variance for PCA
        pcaReady = zeros(size(stackedResps));
        for col = 1:size(stackedResps, 2)
            pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
        end
        
        % Apply PCA and store away variance explained and PCs
        % (outputs I ignore are the principal component coefficients, the variances
        % in the latent space, Hotelling's T-squared statistic for each
        % observation, and the mean of each column)
        [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % I can also go ahead and store max and min xlim, ylim, and zlim
        max_x = max(pcaResult(:, 1)); min_x = min(pcaResult(:, 1));
        max_y = max(pcaResult(:, 2)); min_y = min(pcaResult(:, 2));
        max_z = max(pcaResult(:, 3)); min_z = min(pcaResult(:, 3));
        
        % Separate the PCA scores by group, condition, and subject to make
        % labeling of plots easy (i.e., "unstack")
        plotReady = cell(numGroups, numConds);
        
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                plotReady{grp_idx, cond_idx} = cell(numSubjects(grp_idx), 1);
                for sub_idx = 1:numSubjects(grp_idx)
                    plotReady{grp_idx, cond_idx}{sub_idx} = ...
                        pcaResult(subjLabels == subjectNums{grp_idx}(sub_idx) & ...
                        groupLabels == grp_idx & condLabels == cond_idx, :);
                    
                    % Check size to make sure everything was stored properly
                    if size(plotReady{grp_idx, cond_idx}{sub_idx}, 1) ~= ...
                            numSamples(cond_idx) || ...
                            size(plotReady{grp_idx, cond_idx}{sub_idx}, 2) ...
                            ~= numFeats
                        error('plotReady sizing is incorrect')
                    end
                end
            end
        end
        
        % Branch based on dimensions to plot
        if dim2plot == 2
            % Create figure for separate subplots
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize tight subplot with a subplot per condition arranged
            % horizontally, with each group arranged vertically
            subs = tight_subplot(numGroups, numConds);
            
            % Iterate through groups
            for grp_idx = 1:numGroups
                % Iterate through conditions
                for cond_idx = 1:numConds
                    % Each group and condition will be on its own subplot
                    axes(subs(grp_idx*numConds - numConds + cond_idx));
                    hold on
                    % Iterate through participants and overlay
                    for sub_idx = 1:numSubjects(grp_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                            s = scatter(plotReady{grp_idx, cond_idx}{sub_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 2), ...
                                sizeLabels(t), colorLabels(sum(...
                                numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                    
                    % Set xlim and ylim
                    xlim([min_x, max_x]);
                    ylim([min_y, max_y]);
                    
                    % Label axes and title
                    % Include variance explained by each of the principal components in the
                    % axis labels
                    xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                    if cond_idx == 1    % (I only want y-label on leftmost plot)
                        ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                    end
                    if cond_idx == numConds   % (I only want legend on rightmost plot)
                        legend(subjectIDs{grp_idx}, 'Location', 'eastoutside')
                    end
                    
                    % Title it according to group and condition
                    title([groupNames{grp_idx}, '-', condNames_legend{cond_idx}]);
                    hold off
                end
            end
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            
            % Create figure for overlaid data
            fig = figure(2);
            set(fig, 'Visible', 'on');
            
            % Single subplot per condition
            subs = tight_subplot(1, numConds);
            
            % Iterate through conditions
            for cond_idx = 1:numConds
                % Each condition will be its own subplot
                axes(subs(cond_idx));
                hold on
                for grp_idx = 1:numGroups
                    % Iterate through participants and overlay
                    for sub_idx = 1:numSubjects(grp_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                            s = scatter(plotReady{grp_idx, cond_idx}{sub_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 2), ...
                                sizeLabels(t), colorLabels(...
                                sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
                                :), grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                end
                
                % Set xlim and ylim
                xlim([min_x, max_x]);
                ylim([min_y, max_y]);
                
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                if cond_idx == 1    % (I only want y-label on leftmost plot)
                    ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                end
                if cond_idx == numConds    % (I only want legend on rightmost plot)
                    % Create new subjectIDs_all cell array to use now that we
                    % have all subjects on same plot
                    subjectIDs_all = {};
                    for grp_idx = 1:numGroups
                        subjectIDs_all = [subjectIDs_all; subjectIDs{grp_idx}];
                    end
                    legend(subjectIDs_all, 'Location', 'eastoutside')
                end
                
                % Title it according to condition
                title(condNames_legend{cond_idx});
                
                hold off
                
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Overlaid']);
            
            % Save plot
            saveas(fig, [fileName, ' - Overlaid'], 'fig');
            
            % Close both figures
            close all
            
        else
            % Create figure for separate subplots
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize tight subplot with a subplot per condition arranged
            % horizontally, with each group arranged vertically
            subs = tight_subplot(numGroups, numConds);
            
            % Iterate through groups
            for grp_idx = 1:numGroups
                % Iterate through conditions
                for cond_idx = 1:numConds
                    % Each group and condition will be on its own subplot
                    axes(subs(grp_idx*numConds - numConds + cond_idx));
                    hold on
                    % Iterate through participants and overlay
                    for sub_idx = 1:numSubjects(grp_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                            s = scatter3(plotReady{grp_idx, cond_idx}{sub_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 2), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 3), ...
                                sizeLabels(t), colorLabels(sum(...
                                numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                    
                    % Set xlim and ylim
                    xlim([min_x, max_x]);
                    ylim([min_y, max_y]);
                    zlim([min_z, max_z]);
                    
                    % Label axes and title
                    % Include variance explained by each of the principal components in the
                    % axis labels
                    xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                    ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                    zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                    if cond_idx == numConds   % (I only want legend on rightmost plot)
                        legend(subjectIDs{grp_idx}, 'Location', 'eastoutside')
                    end
                    
                    % Title it according to group and condition
                    title([groupNames{grp_idx}, '-', condNames_legend{cond_idx}]);
                    hold off
                end
            end
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            
            % Create figure for overlaid data
            fig = figure(2);
            set(fig, 'Visible', 'on');
            
            % Single subplot per condition
            subs = tight_subplot(1, numConds);
            
            % Iterate through conditions
            for cond_idx = 1:numConds
                % Each condition will be its own subplot
                axes(subs(cond_idx));
                hold on
                for grp_idx = 1:numGroups
                    % Iterate through participants and overlay
                    for sub_idx = 1:numSubjects(grp_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                            s = scatter3(plotReady{grp_idx, cond_idx}{sub_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 2), ...
                                plotReady{grp_idx, cond_idx}{sub_idx}(t, 3), ...
                                sizeLabels(t), colorLabels(...
                                sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
                                :), grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{sub_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                end
                
                % Set xlim and ylim
                xlim([min_x, max_x]);
                ylim([min_y, max_y]);
                zlim([min_z, max_z]);
                
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                if cond_idx == numConds    % (I only want legend on rightmost plot)
                    % Create new subjectIDs_all cell array to use now that we
                    % have all subjects on same plot
                    subjectIDs_all = {};
                    for grp_idx = 1:numGroups
                        subjectIDs_all = [subjectIDs_all; subjectIDs{grp_idx}];
                    end
                    legend(subjectIDs_all, 'Location', 'eastoutside')
                end
                
                % Title it according to condition
                title(condNames_legend{cond_idx});
                
                hold off
                
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Overlaid']);
            
            % Save plot
            saveas(fig, [fileName, ' - Overlaid'], 'fig');
            
            % Close both figures
            close all
        end
        
        % Create bar graphs for each of the PCs plotted specifying their
        % "loadings" a.k.a. "recipes"
        fig = figure(1);
        set(fig, 'Visible', 'on');
        for dim = 1:dim2plot
            % Create a dim2plot x 1 subplot grid
            subplot(dim2plot, 1, dim)
            
            % Each principal component's loadings is a column of the matrix
            bar(pcaLoads(:, dim), 'red')
            
            % Label axes and title
            xlabel('Feature Index')
            ylabel('Coefficient')
            title(['PC', num2str(dim), ' Recipe'])
        end
        % Set title for entire figure
        sgtitle('Principal Component Loadings')
        
        % Save this figure away
        saveas(fig, [fileName, '_PCloadings'], 'fig');
        
        % Close all plots
        close all

    end
else
    % I'll need to know how many repetitions I have of each condition type
    conds_perType = zeros(numConds, 1);
    for cond_idx = 1:numConds
        % First group and first participant used arbitrarily
        conds_perType(cond_idx) = ...
            size(featSegm_groupSep{1}{1}.(condNames{cond_idx}), 3);
    end
    
    % Separate responses by group and by condition type; within each
    % separation, separate one step further by repetitions of a condition
    % (e.g., neutral 1 vs. neutral 2); within these cells, stack by subject
    resps_groupSep_condSep_repSep = cell(numGroups, numConds);
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            % Initialize
            resps_groupSep_condSep_repSep{grp_idx, cond_idx} = ...
                cell(conds_perType(cond_idx), 1);
            for rep_idx = 1:conds_perType(cond_idx)
                % Initialize
                resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx} = ...
                    zeros(numSamples(cond_idx), numFeats, numSubjects(grp_idx));
            end
        end
    end
    
    % Now fill the nested cell array's arrays
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            for rep_idx = 1:conds_perType(cond_idx)
                for sub = 1:numSubjects(grp_idx)
                    resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}(:, :, sub) = ...
                        featSegm_groupSep{grp_idx}{sub}.(condNames{cond_idx})(:, :, rep_idx);
                end
            end
        end
    end
    
    % To encode time, I am going to create a time vector according to the
    % longest time protocol condition we have and then just use a subset of
    % that for the other conditions
    timeVec = 0:(1/Fs):(max(numSamples)-1)/Fs;
    
    % Encode time in both alpha and size of scatter plot points using
    % alpha and size conversion rule below
    alphaMin = 0.3; sizeMin = 20; timeMin = min(timeVec);
    alphaMax = 1; sizeMax = 70; timeMax = max(timeVec);
    
    % Change in transparency and size per time increment
    alphaSlope = (alphaMax - alphaMin)/(timeMax - timeMin);
    sizeSlope = (sizeMax - sizeMin)/(timeMax - timeMin);
    
    % Store away alpha and size "labels" to encode time
    alphaLabels = alphaSlope*(timeVec - timeMin) + alphaMin;
    sizeLabels = sizeSlope*(timeVec - timeMin) + sizeMin;
    
    % Branch based on meanResp
    if meanResp
        % Set random seed so that colors are always the same
        rng(5); % 0, 1, 2, 3, 4 sucked
        colorLabels = rand(sum(conds_perType), 3);
        
        % To differentiate groups, we will use different markers
        grpLabels = {'o', 'p', 's', '^'};
        
        % Compute means for each condition repetition, within group
        % (3rd dimension of each array is subject)
        meanResps_groupSep_condAvg = cell(numGroups, numConds);
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                % Nest a cell array within to store condition type
                % repetitions
                meanResps_groupSep_condAvg{grp_idx, cond_idx} = ...
                    cell(conds_perType(cond_idx), 1);
                
                % compute means but ignore NaNs since some subjects will
                % not have data for particular condition repetitions
                for rep_idx = 1:conds_perType(cond_idx)
                    meanResps_groupSep_condAvg{grp_idx, cond_idx}{rep_idx} = ...
                        mean(...
                        resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}, ...
                        3, 'omitnan');
                end
            end
        end
        
        % I now need to stack and flatten all of the
        % responses into a 2-D array
        % At the same time, I will need to keep track of groupLabels,
        % condLabels, and repLabels so I can later separate when plotting
        groupLabels = [];
        condLabels = [];
        repLabels = [];
        
        % Initialize the stacked result
        stackedResps = [];
        
        % Loop and stack
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                for rep_idx = 1:conds_perType(cond_idx)
                    % Stack responses
                    stackedResps = [stackedResps; ...
                        meanResps_groupSep_condAvg{grp_idx, cond_idx}{rep_idx}];
                    
                    % Keep track of grouping labels
                    groupLabels = [groupLabels; grp_idx*...
                        ones(size(meanResps_groupSep_condAvg{grp_idx, ...
                        cond_idx}{rep_idx}, 1), 1)];
                    repLabels = [repLabels; rep_idx*...
                        ones(size(meanResps_groupSep_condAvg{grp_idx, ...
                        cond_idx}{rep_idx}, 1), 1)];
                    condLabels = [condLabels; cond_idx*...
                        ones(size(meanResps_groupSep_condAvg{grp_idx, ...
                        cond_idx}{rep_idx}, 1), 1)];
                end
            end
        end
        
        % Check to make sure sizes are correct
        % To do this, I have to calculate how many total samples now that
        % there are repetitions per condition type
        totalSamples = 0;
        for cond_idx = 1:numConds
            totalSamples = totalSamples + ...
                numSamples(cond_idx)*conds_perType(cond_idx);
        end
        
        % Now check size
        if size(stackedResps, 1) ~= numGroups*totalSamples || ...
                length(repLabels) ~=  numGroups*totalSamples || ...
                length(condLabels) ~= numGroups*totalSamples || ...
                length(groupLabels) ~= numGroups*totalSamples
            error('Stacking did not go as planned');
        end
        
        % Now I can make each column unit variance for PCA
        pcaReady = zeros(size(stackedResps));
        for col = 1:size(stackedResps, 2)
            pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
        end
        
        % Apply PCA and store away variance explained and PCs
        % (outputs I ignore are the principal component coefficients, the variances
        % in the latent space, Hotelling's T-squared statistic for each
        % observation, and the mean of each column)
        [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % I can also go ahead and store max and min xlim, ylim, and zlim
        max_x = max(pcaResult(:, 1)); min_x = min(pcaResult(:, 1));
        max_y = max(pcaResult(:, 2)); min_y = min(pcaResult(:, 2));
        max_z = max(pcaResult(:, 3)); min_z = min(pcaResult(:, 3));
        
        
        % Separate the PCA scores by group, condition, and repetition to
        % make labeling of plots easy (i.e., "unstack")
        plotReady = cell(numGroups, numConds);
        
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                % Nested cell arrays to store the repetitions per condition
                plotReady{grp_idx, cond_idx} = cell(conds_perType(cond_idx), 1);
                for rep_idx = 1:conds_perType(cond_idx)
                    plotReady{grp_idx, cond_idx}{rep_idx} = ...
                        pcaResult(repLabels == rep_idx & ...
                        groupLabels == grp_idx & condLabels == cond_idx, :);
                    
                    % Check size to make sure everything was stored properly
                    if size(plotReady{grp_idx, cond_idx}{rep_idx}, 1) ~= ...
                            numSamples(cond_idx) || ...
                            size(plotReady{grp_idx, cond_idx}{rep_idx}, 2) ...
                            ~= numFeats
                        error('plotReady sizing is incorrect')
                    end
                end
            end
        end
        
        if dim2plot == 2
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, numGroups);
            
            % Iterate through groups and condition repetitions and overlay
            for grp_idx = 1:numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:numConds
                    for rep_idx = 1:conds_perType(cond_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{rep_idx}, 1)
                            s = scatter(plotReady{grp_idx, cond_idx}{rep_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{rep_idx}(t, 2), ...
                                sizeLabels(t), ...
                                colorLabels(sum(conds_perType(1:(cond_idx - 1)))...
                                + rep_idx, :), grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{rep_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                
                % Create cell array of strings for legend
                legendNames = cell(sum(conds_perType), 1);
                for cond_idx = 1:numConds
                    for rep_idx = 1:conds_perType(cond_idx)
                        legendNames{sum(conds_perType(1:(cond_idx - 1)))...
                                + rep_idx} = [condNames_legend{cond_idx}, ...
                            num2str(rep_idx)];
                    end
                end
                
                % Set the legend
                legend(legendNames, 'Location', 'best')
                legend('boxoff')
                hold off
                
                % Add title to subplot
                title(groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Because video_lowDtrajectories expects a (# markers) x (#
            % colors) cell array, rather than the nested structure I am
            % using currently, I form a new array that is structured
            % similar to legendNames where the conditions and repetitions
            % are flattened to be one dimension
            videoReady = cell(numGroups, sum(conds_perType));
            for grp_idx = 1:numGroups
                for cond_idx = 1:numConds
                    for rep_idx =1:conds_perType(cond_idx)
                        videoReady{grp_idx, ...
                            sum(conds_perType(1:(cond_idx - 1))) + rep_idx} = ...
                            plotReady{grp_idx, cond_idx}{rep_idx};                  
                    end
                end
            end
            
            % Create video for all data overlaid
            video_lowDtrajectories(videoReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', legendNames, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
            
            % Close all figures
            close all
            
        else
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, numGroups);
            
            % Iterate through groups and conditions and overlay
            for grp_idx = 1:numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:numConds
                    for rep_idx = 1:conds_perType(cond_idx)
                        % Plot one scatter point at a time to leverage time labels
                        for t = 1:size(plotReady{grp_idx, cond_idx}{rep_idx}, 1)
                            s = scatter3(plotReady{grp_idx, cond_idx}{rep_idx}(t, 1), ...
                                plotReady{grp_idx, cond_idx}{rep_idx}(t, 2), ...
                                plotReady{grp_idx, cond_idx}{rep_idx}(t, 3), ...
                                sizeLabels(t), ...
                                colorLabels(sum(conds_perType(1:(cond_idx - 1)))...
                                + rep_idx, :), grpLabels{grp_idx}, 'filled');
                            s.MarkerFaceAlpha = alphaLabels(t);
                            
                            % If this is the final datapoint, then store legend
                            % info
                            if t == size(plotReady{grp_idx, cond_idx}{rep_idx}, 1)
                                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                            else
                                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                
                % Create cell array of strings for legend
                legendNames = cell(sum(conds_perType), 1);
                for cond_idx = 1:numConds
                    for rep_idx = 1:conds_perType(cond_idx)
                        legendNames{sum(conds_perType(1:(cond_idx - 1)))...
                            + rep_idx} = [condNames_legend{cond_idx}, ...
                            num2str(rep_idx)];
                    end
                end
                
                % Set the legend
                legend(legendNames, 'Location', 'best')
                legend('boxoff')
                hold off
                
                % Add title to subplot
                title(groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' - Separate']);
            
            % Save plot
            saveas(fig, [fileName, ' - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Because video_lowDtrajectories expects a (# markers) x (#
            % colors) cell array, rather than the nested structure I am
            % using currently, I form a new array that is structured
            % similar to legendNames where the conditions and repetitions
            % are flattened to be one dimension
            videoReady = cell(numGroups, sum(conds_perType));
            for grp_idx = 1:numGroups
                for cond_idx = 1:numConds
                    for rep_idx =1:conds_perType(cond_idx)
                        videoReady{grp_idx, ...
                            sum(conds_perType(1:(cond_idx - 1))) + rep_idx} = ...
                            plotReady{grp_idx, cond_idx}{rep_idx};                  
                    end
                end
            end
            
            % Create video for all data overlaid
            video_lowDtrajectories(videoReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', legendNames, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
            
            % Close all figures
            close all
        end
        
        % Create bar graphs for each of the PCs plotted specifying their
        % "loadings" a.k.a. "recipes"
        fig = figure(1);
        set(fig, 'Visible', 'on');
        for dim = 1:dim2plot
            % Create a dim2plot x 1 subplot grid
            subplot(dim2plot, 1, dim)
            
            % Each principal component's loadings is a column of the matrix
            bar(pcaLoads(:, dim), 'red')
            
            % Label axes and title
            xlabel('Feature Index')
            ylabel('Coefficient')
            title(['PC', num2str(dim), ' Recipe'])
        end
        % Set title for entire figure
        sgtitle('Principal Component Loadings')
        
        % Save this figure away
        saveas(fig, [fileName, '_PCloadings'], 'fig');
        
        % Close all plots
        close all
        
    else
        % The main challenge to consider here when it comes to plotting the
        % data is that now we will have multiple repetitions per condition
        % type, and all of these repetitions need separate subplots since
        % each subplot will have all the subjects' data overlaid, either
        % from a single group or from the groups being compared. What we
        % will do is separate the condition types with separate figures and
        % have the repetitions have their own subplots. Still, because
        % VNS/sham, for example, has 6 conditions for the DARPA VNS before
        % lunch dataset, we need a max number of subplots per row
        maxPlots_perRow = 3:1:4;   % If you include more than one option, ...
        % ...then the function will choose the option that divides the
        % total repetitions cleanly or leaves the least number of gaps for
        % subplots in the figure. Note that the way this function is
        % written, the options must increment by one and be contiguous
        % whole numbers
        
        
        % Create color map based on total number of subjects
        rng(0);     % Set random seed so that colors are always the same
        colorLabels = rand(sum(numSubjects), 3);
        
        % To differentiate groups, we will use different markers
        grpLabels = {'o', 'p', 's', '^'};
        
        % Create nested cell array of subject IDs based on groups and subject
        % numbers
        subjectIDs = cell(numGroups, 1);
        for grp_idx = 1:numGroups
            subjectIDs{grp_idx} = cell(numSubjects(grp_idx), 1);
        end
        for grp_idx = 1:numGroups
            for sub = 1:numel(subjectIDs{grp_idx})
                subjectIDs{grp_idx}{sub} = num2str(subjectNums{grp_idx}(sub));
            end
        end
        
        % I now need to stack and flatten all of the
        % responses into a 2-D array
        % At the same time, I will need to keep track of groupLabels,
        % condLabels, repLabels, and subjLabels so I can later separate 
        % when plotting
        groupLabels = [];
        condLabels = [];
        repLabels = [];
        subjLabels = [];
        
        % Initialize the stacked result
        stackedResps = [];
        
        % Loop and stack
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                for rep_idx = 1:conds_perType(cond_idx)
                    for sub_idx = 1:numSubjects(grp_idx)
                        % Stack responses, unless they are NaN
                        if ~isnan(sum(resps_groupSep_condSep_repSep{grp_idx, ...
                                cond_idx}{rep_idx}(:, :, sub_idx), 'all'))
                            stackedResps = [stackedResps; ...
                                resps_groupSep_condSep_repSep{grp_idx, ...
                                cond_idx}{rep_idx}(:, :, sub_idx)];
                            
                            % Keep track of grouping labels
                            groupLabels = [groupLabels; grp_idx*...
                                ones(size(...
                                resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}...
                                (:, :, sub_idx), 1), 1)];
                            condLabels = [condLabels; cond_idx*...
                                ones(size(...
                                resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}...
                                (:, :, sub_idx), 1), 1)];
                            repLabels = [repLabels; rep_idx*...
                                ones(size(...
                                resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}...
                                (:, :, sub_idx), 1), 1)];
                            subjLabels = [subjLabels; subjectNums{grp_idx}(sub_idx)*...
                                ones(size(...
                                resps_groupSep_condSep_repSep{grp_idx, cond_idx}{rep_idx}...
                                (:, :, sub_idx), 1), 1)];
                        end
                    end
                end
            end
        end
        
        % Now I can make each column unit variance for PCA
        pcaReady = zeros(size(stackedResps));
        for col = 1:size(stackedResps, 2)
            pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
        end
        
        % Apply PCA and store away variance explained and PCs
        % (outputs I ignore are the principal component coefficients, the variances
        % in the latent space, Hotelling's T-squared statistic for each
        % observation, and the mean of each column)
        [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % I can also go ahead and store max and min xlim, ylim, and zlim
        max_x = max(pcaResult(:, 1)); min_x = min(pcaResult(:, 1));
        max_y = max(pcaResult(:, 2)); min_y = min(pcaResult(:, 2));
        max_z = max(pcaResult(:, 3)); min_z = min(pcaResult(:, 3));
        
        % Separate the PCA scores by group, condition, repetition, and
        % subject to make labeling of plots easy (i.e., "unstack")
        plotReady = cell(numGroups, numConds);
        
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                % Nested cell arrays to store the repetitions per condition
                plotReady{grp_idx, cond_idx} = cell(conds_perType(cond_idx), 1);
                for rep_idx = 1:conds_perType(cond_idx)
                    % Nested cell arrays to store subjects' data separately
                    plotReady{grp_idx, cond_idx}{rep_idx} = cell(numSubjects(grp_idx), 1);
                    for sub_idx = 1:numSubjects(grp_idx)
                        % If data does not exist, then the cell will remain
                        % empty, and I can find this using isempty()
                        plotReady{grp_idx, cond_idx}{rep_idx}{sub_idx} = ...
                            pcaResult(subjLabels == subjectNums{grp_idx}(sub_idx) & ...
                            groupLabels == grp_idx & condLabels == cond_idx ...
                            & repLabels == rep_idx, :);
                        
                        % Check size to make sure everything was stored properly
                        if ~isempty(plotReady{grp_idx, cond_idx}{rep_idx}{sub_idx}) ...
                                && (size(plotReady{grp_idx, cond_idx}{...
                                rep_idx}{sub_idx}, 1) ~= numSamples(cond_idx) || ...
                                size(plotReady{grp_idx, cond_idx}{rep_idx}{sub_idx}, 2) ...
                                ~= numFeats)
                            error('plotReady sizing is incorrect')
                        end
                    end
                end
            end
        end
        
        % Branch based on dimension to plot
        if dim2plot == 2
            % We are going to create a figure per condition type because
            % there are too many total condition repetitions to plot it all
            % in one figure with subjects' data overlaid
            for cond_idx = 1:numConds
                % Create figure for separate subplots for separate groups
                fig = figure(1);
                set(fig, 'Visible', 'on');
                
                % Figure out how many rows we will need
                if conds_perType(cond_idx) <= maxPlots_perRow(end)
                    subs_perRow = conds_perType(cond_idx);
                else
                    perfectFound = false;   % Divides evenly
                    subs_perRow = maxPlots_perRow(end); % What number to use
                    bestDiff = inf;   % Helps find which number to use
                    % Search
                    for optn = flip(maxPlots_perRow)
                        % Only do the following if we haven't already found
                        % the perfect answer
                        if ~perfectFound
                            % If this divides perfectly, set the
                            % perfectFound boolean to true to effectivey
                            % stop the search
                            if mod(conds_perType(cond_idx), optn) == 0
                                perfectFound = true;
                                subs_perRow = optn;
                                bestDiff = 0;
                            else
                                % Search for the smallest difference
                                % between the subplots per row option and
                                % the remainder; that will help us minimize
                                % the gaps in the subplot grid
                                if optn - mod(conds_perType(cond_idx), optn) < ...
                                        bestDiff
                                    bestDiff = ...
                                        optn - mod(conds_perType(cond_idx), optn);
                                    subs_perRow = optn;
                                end
                            end
                        end
                    end
                end
                
                % Initialize subplot grid
                numRows_perGrp = ceil(conds_perType(cond_idx)/subs_perRow);
                subs = tight_subplot(numGroups*numRows_perGrp, subs_perRow);
                
                plot_idx = 1;    % Counter to help plot properly
                % Iterate through groups
                for grp_idx = 1:numGroups
                    % Iterate through repetitions of the condition
                    for rep_idx = 1:conds_perType(cond_idx)
                        % I will only add the subject IDs for those
                        % subjects with data
                        legendNames = {};
                        
                        % Each group and repetition will be on its own
                        % subplot for this figure
                        axes(subs(plot_idx));
                        hold on
                        % Iterate through participants and overlay
                        for sub_idx = 1:numSubjects(grp_idx)
                            % Plot one scatter point at a time to leverage time labels
                            if ~isempty(plotReady{grp_idx, cond_idx}{...
                                    rep_idx}{sub_idx})
                                % Add this subject ID to legend
                                legendNames = [legendNames; ...
                                    subjectIDs{grp_idx}{sub_idx}];
                                for t = 1:size(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}, 1)
                                    s = scatter(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 1), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 2), ...
                                        sizeLabels(t), colorLabels(sum(...
                                        numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                        grpLabels{grp_idx}, 'filled');
                                    s.MarkerFaceAlpha = alphaLabels(t);
                                    
                                    % If this is the final datapoint, then store legend
                                    % info
                                    if t == size(plotReady{grp_idx, cond_idx}{...
                                            rep_idx}{sub_idx}, 1)
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                                    else
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                                    end
                                end
                            end
                        end
                        
                        % Set xlim and ylim
                        xlim([min_x, max_x]);
                        ylim([min_y, max_y]);
                        
                        % Label axes and title
                        % Include variance explained by each of the principal components in the
                        % axis labels
                        xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                        ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                        
                        % Set legend only based on participants with data
                        legend(legendNames)
                        
                        % Title it according to group and condition
                        % repetition
                        title([groupNames{grp_idx}, '-', condNames_legend{cond_idx}, ...
                            '-', num2str(rep_idx)]);
                        hold off
                        
                        % Increment counter
                        if rep_idx == conds_perType(cond_idx) && ...
                                mod(plot_idx, subs_perRow) ~= 0
                            % We enter this section of code only if we are
                            % at the end of condition repetitions for a
                            % group but we are not at the rightmost edge of
                            % the subplot grid yet. That means, we need to
                            % skip over a gap to go onto the next row for
                            % our next plots
                            plot_idx = plot_idx + (subs_perRow - ...
                                mod(plot_idx, subs_perRow)) + 1;
                        else
                            plot_idx = plot_idx + 1;
                        end
                    end
                end
                
                % Add title to entire figure
                sgtitle([plotTitle, ' - Separate - Condition Type: ', ...
                    condNames_legend{cond_idx}]);
                
                % Save plot
                saveas(fig, [fileName, ' - Separate_', ...
                    condNames{cond_idx}], 'fig');
                
                
                
                % Create figure for groups overlaid
                fig = figure(2);
                set(fig, 'Visible', 'on');
                
                subs = tight_subplot(numRows_perGrp, subs_perRow);
                plot_idx = 1;    % Counter to help plot properly
                % Iterate through condition repetitions
                for rep_idx = 1:conds_perType(cond_idx)
                    % I will only add the subject IDs for those
                    % subjects with data
                    legendNames = {};
                    
                    axes(subs(plot_idx));
                    hold on
                    % Iterate through groups
                    for grp_idx = 1:numGroups
                        % Iterate through participants and overlay
                        for sub_idx = 1:numSubjects(grp_idx)
                            % Plot one scatter point at a time to leverage time labels
                            if ~isempty(plotReady{grp_idx, cond_idx}{...
                                    rep_idx}{sub_idx})
                                % Add this subject ID to legend
                                legendNames = [legendNames; ...
                                    subjectIDs{grp_idx}{sub_idx}];
                                for t = 1:size(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}, 1)
                                    s = scatter(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 1), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 2), ...
                                        sizeLabels(t), colorLabels(sum(...
                                        numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                        grpLabels{grp_idx}, 'filled');
                                    s.MarkerFaceAlpha = alphaLabels(t);
                                    
                                    % If this is the final datapoint, then store legend
                                    % info
                                    if t == size(plotReady{grp_idx, cond_idx}{...
                                            rep_idx}{sub_idx}, 1)
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                                    else
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                                    end
                                end
                            end
                        end
                    end
                    
                    % Set xlim and ylim
                    xlim([min_x, max_x]);
                    ylim([min_y, max_y]);
                    
                    % Label axes and title
                    % Include variance explained by each of the principal components in the
                    % axis labels
                    xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                    ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                    
                    % Set legend only based on participants with data
                    legend(legendNames)
                    
                    % Title it according to condition repetition
                    title([condNames_legend{cond_idx}, '-', num2str(rep_idx)]);
                    hold off
                    
                    % Increment counter
                    plot_idx = plot_idx + 1;
                end
                
                % Add title to entire figure
                sgtitle([plotTitle, ' - Overlaid - Condition Type: ', ...
                    condNames_legend{cond_idx}]);
                
                % Save plot
                saveas(fig, [fileName, ' - Overlaid_', ...
                    condNames{cond_idx}], 'fig');
                
                % Close both figures
                close all
            end
            
        else
            % We are going to create a figure per condition type because
            % there are too many total condition repetitions to plot it all
            % in one figure with subjects' data overlaid
            for cond_idx = 1:numConds
                % Create figure for separate subplots for separate groups
                fig = figure(1);
                set(fig, 'Visible', 'on');
                
                % Figure out how many rows we will need
                if conds_perType(cond_idx) <= maxPlots_perRow(end)
                    subs_perRow = conds_perType(cond_idx);
                else
                    perfectFound = false;   % Divides evenly
                    subs_perRow = maxPlots_perRow(end); % What number to use
                    bestDiff = inf;   % Helps find which number to use
                    % Search
                    for optn = flip(maxPlots_perRow)
                        % Only do the following if we haven't already found
                        % the perfect answer
                        if ~perfectFound
                            % If this divides perfectly, set the
                            % perfectFound boolean to true to effectivey
                            % stop the search
                            if mod(conds_perType(cond_idx), optn) == 0
                                perfectFound = true;
                                subs_perRow = optn;
                                bestDiff = 0;
                            else
                                % Search for the smallest difference
                                % between the subplots per row option and
                                % the remainder; that will help us minimize
                                % the gaps in the subplot grid
                                if optn - mod(conds_perType(cond_idx), optn) < ...
                                        bestDiff
                                    bestDiff = ...
                                        optn - mod(conds_perType(cond_idx), optn);
                                    subs_perRow = optn;
                                end
                            end
                        end
                    end
                end
                
                % Initialize subplot grid
                numRows_perGrp = ceil(conds_perType(cond_idx)/subs_perRow);
                subs = tight_subplot(numGroups*numRows_perGrp, subs_perRow);
                
                plot_idx = 1;    % Counter to help plot properly
                % Iterate through groups
                for grp_idx = 1:numGroups
                    % Iterate through repetitions of the condition
                    for rep_idx = 1:conds_perType(cond_idx)
                        % I will only add the subject IDs for those
                        % subjects with data
                        legendNames = {};
                        
                        % Each group and repetition will be on its own
                        % subplot for this figure
                        axes(subs(plot_idx));
                        hold on
                        % Iterate through participants and overlay
                        for sub_idx = 1:numSubjects(grp_idx)
                            % Plot one scatter point at a time to leverage time labels
                            if ~isempty(plotReady{grp_idx, cond_idx}{...
                                    rep_idx}{sub_idx})
                                % Add this subject ID to legend
                                legendNames = [legendNames; ...
                                    subjectIDs{grp_idx}{sub_idx}];
                                for t = 1:size(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}, 1)
                                    s = scatter3(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 1), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 2), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 3), ...
                                        sizeLabels(t), colorLabels(sum(...
                                        numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                        grpLabels{grp_idx}, 'filled');
                                    s.MarkerFaceAlpha = alphaLabels(t);
                                    
                                    % If this is the final datapoint, then store legend
                                    % info
                                    if t == size(plotReady{grp_idx, cond_idx}{...
                                            rep_idx}{sub_idx}, 1)
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                                    else
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                                    end
                                end
                            end
                        end
                        
                        % Set xlim and ylim
                        xlim([min_x, max_x]);
                        ylim([min_y, max_y]);
                        zlim([min_z, max_z]);
                        
                        % Label axes and title
                        % Include variance explained by each of the principal components in the
                        % axis labels
                        xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                        ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                        zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                        
                        % Set legend only based on participants with data
                        legend(legendNames)
                        
                        % Title it according to group and condition
                        % repetition
                        title([groupNames{grp_idx}, '-', condNames_legend{cond_idx}, ...
                            '-', num2str(rep_idx)]);
                        hold off
                        
                        % Increment counter
                        if rep_idx == conds_perType(cond_idx) && ...
                                mod(plot_idx, subs_perRow) ~= 0
                            % We enter this section of code only if we are
                            % at the end of condition repetitions for a
                            % group but we are not at the rightmost edge of
                            % the subplot grid yet. That means, we need to
                            % skip over a gap to go onto the next row for
                            % our next plots
                            plot_idx = plot_idx + (subs_perRow - ...
                                mod(plot_idx, subs_perRow)) + 1;
                        else
                            plot_idx = plot_idx + 1;
                        end
                    end
                end
                
                % Add title to entire figure
                sgtitle([plotTitle, ' - Separate - Condition Type: ', ...
                    condNames_legend{cond_idx}]);
                
                % Save plot
                saveas(fig, [fileName, ' - Separate_', ...
                    condNames{cond_idx}], 'fig');
                
                
                
                % Create figure for groups overlaid
                fig = figure(2);
                set(fig, 'Visible', 'on');
                
                subs = tight_subplot(numRows_perGrp, subs_perRow);
                plot_idx = 1;    % Counter to help plot properly
                % Iterate through condition repetitions
                for rep_idx = 1:conds_perType(cond_idx)
                    % I will only add the subject IDs for those
                    % subjects with data
                    legendNames = {};
                    
                    axes(subs(plot_idx));
                    hold on
                    % Iterate through groups
                    for grp_idx = 1:numGroups
                        % Iterate through participants and overlay
                        for sub_idx = 1:numSubjects(grp_idx)
                            % Plot one scatter point at a time to leverage time labels
                            if ~isempty(plotReady{grp_idx, cond_idx}{...
                                    rep_idx}{sub_idx})
                                % Add this subject ID to legend
                                legendNames = [legendNames; ...
                                    subjectIDs{grp_idx}{sub_idx}];
                                for t = 1:size(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}, 1)
                                    s = scatter3(plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 1), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 2), ...
                                        plotReady{grp_idx, cond_idx}{...
                                        rep_idx}{sub_idx}(t, 3), ...
                                        sizeLabels(t), colorLabels(sum(...
                                        numSubjects(1:(grp_idx-1))) + sub_idx, :), ...
                                        grpLabels{grp_idx}, 'filled');
                                    s.MarkerFaceAlpha = alphaLabels(t);
                                    
                                    % If this is the final datapoint, then store legend
                                    % info
                                    if t == size(plotReady{grp_idx, cond_idx}{...
                                            rep_idx}{sub_idx}, 1)
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                                    else
                                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                                    end
                                end
                            end
                        end
                    end
                    
                    % Set xlim and ylim
                    xlim([min_x, max_x]);
                    ylim([min_y, max_y]);
                    zlim([min_z, max_z]);
                    
                    % Label axes and title
                    % Include variance explained by each of the principal components in the
                    % axis labels
                    xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                    ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                    zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
                    
                    % Set legend only based on participants with data
                    legend(legendNames)
                    
                    % Title it according to condition repetition
                    title([condNames_legend{cond_idx}, '-', num2str(rep_idx)]);
                    hold off
                    
                    % Increment counter
                    plot_idx = plot_idx + 1;
                end
                
                % Add title to entire figure
                sgtitle([plotTitle, ' - Overlaid - Condition Type: ', ...
                    condNames_legend{cond_idx}]);
                
                % Save plot
                saveas(fig, [fileName, ' - Overlaid_', ...
                    condNames{cond_idx}], 'fig');
                
                % Close both figures
                close all
            end
            
        end
        
        
        % Create bar graphs for each of the PCs plotted specifying their
        % "loadings" a.k.a. "recipes"
        fig = figure(1);
        set(fig, 'Visible', 'on');
        for dim = 1:dim2plot
            % Create a dim2plot x 1 subplot grid
            subplot(dim2plot, 1, dim)
            
            % Each principal component's loadings is a column of the matrix
            bar(pcaLoads(:, dim), 'red')
            
            % Label axes and title
            xlabel('Feature Index')
            ylabel('Coefficient')
            title(['PC', num2str(dim), ' Recipe'])
        end
        % Set title for entire figure
        sgtitle('Principal Component Loadings')
        
        % Save this figure away
        saveas(fig, [fileName, '_PCloadings'], 'fig');
        
        % Close all plots
        close all
        
    end
    
end


end