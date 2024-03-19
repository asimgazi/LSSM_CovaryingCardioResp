function plot_realResponses_groupSpec(groupName, saveDir, featSegm_group, ...
    Fs, varargin)
% This function takes segmented responses and plots and saves them
% according to inputs provided

% --------- Inputs --------------- %
% featPer: number of feature subplots per figure
% saveDir: string specifying the directory where plots should be saved
% featSegm_group: Nx1 cell array, where N is number of subjects in group.
%                 Each cell contains a struct that contains a field
%                 for each protocol condition, where within each field, an
%                 array exists that is TxMxC, where T is the number of time
%                 points in that condition segment, M is the number of
%                 features, and C is the number of condition repetitions
%                 for that condition
% Fs: sampling frequency; used to construct time vectors to plot against
% (optional)
% 'plotOption': flag that specifies to look for the option for plotting in
%               next function input. Options to choose from:
%               'condSep': each protocol condition separately (i.e.,
%               neutral 1 is plotted separate from neutral 2, etc.). If a
%               number is included as the next argument after 'condSep',
%               then that number is assumed to be featPer
%               'condAvg': if this option is selected, then the next input
%               should be the number of features per plot since conditions
%               are overlaid on same subplot
%               'condAvg_train': if this option is selected, then the next
%               two inputs should be number of features and the train-test
%               split condition code (e.g., 6.5 for split between 6 and 7)
% 'names': flag that specifies next input as cell array of feature names.
%          If no names are provided, the responses are labeled by index
% 'subjOverlay' - if this flag is provided, the next varargin cell needs to
% be an vector of subject numbers to use for a legend. What this flag
% entails is it will set meanSEM = false. This will overlay all subjects'
% responses for a particular condition on top of each other while
% separating conditions into separate plots. If meanSEM = true, mean +- SEM
% will be computed across participants for each condition and then overlaid
% on the same plot for easy comparison, so the legend will just be based on
% condition as it was for subject-specific plots

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'condSep')
        condAvg = false;
        trainOnly = false;
        % If the next argument is a number, then that means the user wants
        % a specific featPer
        if isnumeric(varargin{arg + 1})
            featPer = varargin{arg + 1};
        end
    elseif strcmp(varargin{arg}, 'condAvg')
        condAvg = true;
        trainOnly = false;
        featPer = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'condAvg_train')
        condAvg = true;
        trainOnly = true;
        featPer = varargin{arg + 1};
        trTest_split = varargin{arg + 2};
    elseif strcmp(varargin{arg}, 'names')
        featNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanSEM = false;    % We will plot responses themselves
        subjectNums = varargin{arg + 1};    % Subject numbers for legend
    end
end

% Set defaults
if ~exist('condAvg', 'var')
    condAvg = true;
end
if ~exist('trainOnly', 'var')
    trainOnly = true;
end
if ~exist('featPer', 'var')
    featPer = 6;
end
if ~exist('trTest_split', 'var')
    trTest_split = 6.5;
end
if ~exist('featNames', 'var')
    % If they don't provide feature names, I'll just number them
    featNames = {};
    for m = 1:size(featSegm_group{1}.n, 2)
        featNames = [featNames; {['Feat. ', num2str(m)]}];
    end
end
if ~exist('meanSEM', 'var')
    meanSEM = true; % read documentation above for what this bool dictates
end

% Just so we have them as easy to understand variables
numSubjects = numel(featSegm_group);
numFeats = size(featSegm_group{1}.n, 2);
condNames = fieldnames(featSegm_group{1});
numConds = length(condNames);

if condAvg
    % Find out how many features we're dealing with and figures we'll need
    numFigs = ceil(numFeats/featPer);
    
    % Initialize a data structure to hold condition-averaged responses for
    % each subject
    feat_condAvg_group = cell(numSubjects, 1);
    for sub = 1:numSubjects
        feat_condAvg_group{sub} = struct();
        for cond = 1:numConds
            % Initialized to NaNs as default in case no data are available
            % to replace (will ignore NaN responses during plotting)
            feat_condAvg_group{sub}.(condNames{cond}) = ...
                nan*ones(size(featSegm_group{sub}.(condNames{cond}), 1), ...
                numFeats);
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
        for sub = 1:numSubjects
            for cond = 1:numConds
                % Average along the 3rd dimension to average over repetitions
                % Only average the training data
                feat_condAvg_group{sub}.(condNames{cond}) = ...
                    mean(featSegm_group{sub}.(condNames{cond})(:, :, ...
                    1:trainNums(cond)), 3, 'omitnan');
            end
        end
    else
        % In this case, we can just average over whatever we have since no
        % train/test considerations
        for sub = 1:numSubjects
            for cond = 1:numConds
                % Average along the 3rd dimension to average over repetitions
                feat_condAvg_group{sub}.(condNames{cond}) = ...
                    mean(featSegm_group{sub}.(condNames{cond}), 3, 'omitnan');
            end
        end
    end
    
    
    % Plotting code branches off based on meanSEM boolean variable
    if meanSEM
        % Create basic color labels (blue, red, green) for protocol conditions
        colorLabels = (1/255)*[221, 85, 255; ...
            255, 85, 85; ...
            95, 211, 188];
        
        % Since I will need to be working with mean and SEM, I need to go ahead
        % and compute that information for each of the conditions. Let's
        % initialize a cell array, one cell per condition, to make going about
        % this easier (agg is short for aggregate)
        agg_responses = cell(numConds, 1);
        for cond = 1:numConds
            for sub = 1:numSubjects
                % Make sure it's not NaN before stacking in 3rd dimension
                if ~isnan(sum(feat_condAvg_group{sub}.(condNames{cond}), 'all'))
                    agg_responses{cond} = ...
                        cat(3, agg_responses{cond}, ...
                        feat_condAvg_group{sub}.(condNames{cond}));
                end
            end
        end
        
        % Compute mean and SEM for each feature for each time point
        mean_responses = cell(numConds, 1);
        SEM_responses = cell(numConds, 1);
        for cond_idx = 1:numConds
            mean_responses{cond_idx} = mean(agg_responses{cond_idx}, 3);
            SEM_responses{cond_idx} = (1/sqrt(numSubjects))*...
                std(agg_responses{cond_idx}, 0, 3);
        end
        
        % Now I can start creating figures
        feat_idx = 1;           % To help us iterate through features
        for fig_idx = 1:numFigs
            % Initialize figure
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize a tight subplot
            subs = tight_subplot(featPer, 1);
            
            % Iterate through subplots
            for subP_idx = 1:featPer
                % Check to make sure we still have features left to plot
                if feat_idx <= numFeats
                    % Plot in subplot
                    axes(subs(subP_idx));
                    
                    % Iterate through conditions and overlay
                    hold on
                    for cond_idx = 1:numConds
                        % Since conditions have different time lengths,
                        % create a timeVec for each condition
                        timeVec = 0:(1/Fs):...
                            (size(mean_responses{cond_idx}, 1) - 1)/Fs;
                        shadedBounds_linePlot(...
                            mean_responses{cond_idx}(:, feat_idx), ...
                            mean_responses{cond_idx}(:, feat_idx) + ...
                            SEM_responses{cond_idx}(:, feat_idx), ...
                            mean_responses{cond_idx}(:, feat_idx) - ...
                            SEM_responses{cond_idx}(:, feat_idx), ...
                            timeVec, colorLabels(cond_idx, :));
                    end
                    % Label y axis
                    ylabel(['\Delta ', featNames{feat_idx}]);
                    
                    % I only want xlabels on the bottom plots
                    if subP_idx == featPer || feat_idx == numFeats
                        xlabel('Time (s)');
                    end
                    
                    % Add legend
                    legend(condNames);
                    hold off
                end
                
                % Increment feature index
                feat_idx = feat_idx + 1;
            end
            
            % Add title to entire subplot grid
            sgtitle([groupName, ' Mean \pm SEM Figure ', num2str(fig_idx), ...
                ' of ', num2str(numFigs)]);
            
            % Save plot
            saveas(fig, [saveDir, groupName, '_Mean_SEM_', ...
                num2str(fig_idx)], 'fig');
            close all
        end
    else
        % Create color map based on number of subjects
        rng(0);     % Set random seed so that colors are always the same
        colorLabels = rand(numSubjects, 3);
        
        % Create cell array of subject IDs based on subject numbers
        subjectIDs = cell(length(subjectNums), 1);
        for sub = 1:numel(subjectIDs)
            subjectIDs{sub} = num2str(subjectNums(sub));
        end
        
        % Start creating figures
        feat_idx = 1;   % To help us iterate through features
        for fig_idx = 1:numFigs
            % Initialize figure
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize a tight subplot grid
            subs = tight_subplot(featPer, numConds);
            
            % Loop through features in this fig. (vertically through subplots)
            for row_idx = 1:featPer
                % Check to make sure we still have features left to plot
                if feat_idx <= numFeats
                    % So I can make y-axis limits all the same for a particular
                    % feature, I need to find out the maximum and minimum
                    % feature values across all participants
                    max_perSubj_perCond = zeros(numSubjects, numConds);
                    min_perSubj_perCond = zeros(numSubjects, numConds);
                    max_perSubj = zeros(numSubjects, 1);
                    min_perSubj = zeros(numSubjects, 1);
                    for sub = 1:numSubjects
                        for cond = 1:numConds
                            max_perSubj_perCond(sub, cond) = max(feat_condAvg_group{...
                                sub}.(condNames{cond})(:, feat_idx));
                            min_perSubj_perCond(sub, cond) = min(feat_condAvg_group{...
                                sub}.(condNames{cond})(:, feat_idx));
                        end
                        max_perSubj = max(max_perSubj_perCond, [], 2);
                        min_perSubj = min(min_perSubj_perCond, [], 2);
                    end
                    
                    % Figure out ylim min and max
                    max_yaxis = max(max_perSubj);
                    min_yaxis = min(min_perSubj);
                    
                    % Loop through conditions (horizontally through subplots)
                    for cond_idx = 1:numConds
                        % Plot in appropriate subplot
                        axes(subs(row_idx*numConds - numConds + cond_idx));
                        
                        % Iterate through participants and overlay
                        hold on
                        for sub = 1:numSubjects
                            % Since conditions have different time lengths,
                            % just create time vectors separately
                            timeVec = 0:(1/Fs):...
                                (size(feat_condAvg_group{sub}.(condNames{cond_idx}), 1)...
                                - 1)/Fs;
                            plot(timeVec, ...
                                feat_condAvg_group{sub}.(condNames{cond_idx})(:, feat_idx), ...
                                'Color', colorLabels(sub, :));
                        end
                        
                        % Set ylim
                        ylim([min_yaxis, max_yaxis])
                        
                        % I only want a ylabel on the leftmost plots
                        if cond_idx == 1
                            ylabel(featNames{feat_idx});
                        end
                        % I only want xlabels on the bottom plots
                        if row_idx == featPer || feat_idx == numFeats
                            xlabel('Time (s)');
                        end
                        % I only want subplot titles on the top plots
                        if row_idx == 1
                            title(condNames{cond_idx});
                        end
                        
                        % I want to create one legend for all, so I will just
                        % arbitrarily assign one of the subplots on the right
                        % to have a legend on the outside bottom right
                        if row_idx == 1 && cond_idx == numConds
                            legend(subjectIDs, 'Location', 'southeastoutside')
                        end
                        
                        hold off
                    end
                    
                    % Increment feature index
                    feat_idx = feat_idx + 1;
                end
            end
            
            % Add title to entire subplot grid
            sgtitle([groupName, ' Subject Overlay Figure ', num2str(fig_idx), ...
                ' of ', num2str(numFigs)]);
            
            % Save plot
            saveas(fig, [saveDir, groupName, '_subjOverlay_', ...
                num2str(fig_idx)], 'fig');
            close all
        end
        
    end
    
else
    % Same two options need to be handled for meanSEM; otherwise, this
    % should be less code since you don't need to average over conditions
    % of the same type anymore. Here, we are treating each condition
    % separately and comparing the protocol conditions themselves rather
    % than aggregating based on what type of condition. Essentially,
    % instead of thinking of neutral vs. trauma, we're thinking condition 1
    % (neutral) vs. condition 2 (neutral) vs. condition 3 (trauma) vs. etc.
    
    % The part that will actually take some figuring out is how do I plan
    % to organize the figures. Because before, we knew that there would
    % only be so many conditions on a single plot, but now, I will have
    % around 10 conditions just for before lunch DARPA VNS data. I can't
    % just overlay them all on top of each other for meanSEM == true plots
    % and I can't just place them next to each other for meanSEM == false
    % plots because 10 subplots horizontally would be ridiculous...
    
    % For meanSEM == false, I am thinking what I could do is create a 
    % figure per feature. Each condition type has an assigned set of rows, 
    % and those assigned rows get filled with each protocol condition of 
    % that type. For the DARPA VNS dataset, there's a max of 14, so if I do
    % two columns, that'd make 7 rows
    
    % For meanSEM == true, I am actually going to overlay protocol
    % conditions of the same type on top of each other, so I will have
    % multiple features per figure, as before, with featPer dictating how
    % many. But then I will have columns == number of condition types,
    % where conditions of the same type will be overlaid on the same plot
    % using different colors. This will help me compare the conditions of
    % the same type of each other (e.g., to answer the question, is the
    % first neutral response same as the other neutral responses or more
    % similar to traumatic stress?) while also being able to compare the
    % general trends across conditions between the condition types
    
    % Plotting code branches off based on meanSEM boolean variable
    if meanSEM
%         % Create color map based on maximum number of conditions for a
%         % condition type
        conds_perType = zeros(numConds, 1);
        for cond_idx = 1:numConds
            % First participant used arbitrarily
            conds_perType(cond_idx) = size(featSegm_group{1}.(condNames{cond_idx}), 3);
        end
%         maxConds_perType = max(conds_perType);
%         % Set random seed so that colors are always the same
%         rng(5); % 0, 1, 2, 3, 4 sucked
%         colorLabels = rand(maxConds_perType, 3);
        
        % Set colors according to condition and repetition
        % Stim. colors
        colorLabels{1} = (1/255)*...
            [141, 95, 211; ...
            212, 42, 255; ...
            255, 85, 221; ...
            68, 33, 120; ...
            222, 135, 205; ...
            85, 0, 212];
        
        % Trauma colors
        colorLabels{2} = (1/255)*...
            [255, 42, 127; ...
            255, 85, 85; ...
            211, 95, 95; ...
            211, 95, 141];
        
        % Neutral colors
        colorLabels{3} = (1/255)*...
            [85, 212, 0; ...
            135, 222, 135; ...
            44, 160, 90; ...
            95, 211, 188];
        
        % SEt marker labels according to condition repetition
        markerLabels = {'o', 'p', 's', '^', 'd', 'x'};
        
        % Since I will need to be working with mean and SEM, I need to go 
        % ahead and compute that information for each of the conditions.
        % To aggregate, I have to keep in mind that not all participants
        % will have data for each of the protocol conditions. To keep
        % things clean and compute SEMs in a kosher way (will need to
        % divide by number of subjects used when computing standard
        % deviation), I will actually separate all protocol conditions out,
        % aggregate, compute means and SEMs, and then reorganize back into
        % separated by condition type using a struct so then it's easy to
        % know what needs to be plotted on top of each other
        agg_responses = cell(sum(conds_perType), 1);
        cond = 1; cond_count = 1;
        for iter = 1:sum(conds_perType)
            for sub = 1:numSubjects
                % Make sure it's not NaN before stacking in 3rd dimension
                if ~isnan(sum(featSegm_group{sub}.(...
                        condNames{cond})(:, :, cond_count), 'all'))
                    agg_responses{iter} = ...
                        cat(3, agg_responses{iter}, ...
                        featSegm_group{sub}.(condNames{cond})(:, :, cond_count));
                end
            end
            
            % Increment/reset other counters
            if cond_count == conds_perType(cond)
                cond_count = 1;
                cond = cond + 1;
            else
                cond_count = cond_count + 1;
            end
            
        end
        
        % Compute mean and SEM for each feature for each time point
        mean_responses = cell(sum(conds_perType), 1);
        SEM_responses = cell(sum(conds_perType), 1);
        for iter = 1:sum(conds_perType)
            mean_responses{iter} = mean(agg_responses{iter}, 3);
            SEM_responses{iter} = (1/sqrt(size(agg_responses{iter}, 3)))*...
                std(agg_responses{iter}, 0, 3);
        end
        
        % Reorganize into struct with separate fields for the types of
        % conditions
        plotReady_mean = struct(); plotReady_SEM = struct();
        iter = 1;
        for cond_idx = 1:numConds
            for cond_count = 1:conds_perType(cond_idx)
                plotReady_mean.(condNames{cond_idx})(:, :, cond_count) = ...
                    mean_responses{iter};
                plotReady_SEM.(condNames{cond_idx})(:, :, cond_count) = ...
                    SEM_responses{iter};
                
                % Increment counter
                iter = iter + 1;
            end
        end
        
        % Find out how many features we're dealing with and figures we'll need
        numFigs = ceil(numFeats/featPer);
        
        % Now I can start creating figures
        feat_idx = 1;   % To help us iterate through features
        for fig_idx = 1:numFigs
            % Initialize figure
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize a tight subplot grid
            subs = tight_subplot(featPer, numConds);
            
            % Loop through features in this fig. (vertically through subplots)
            for row_idx = 1:featPer
                % Check to make sure we still have features left to plot
                if feat_idx <= numFeats
                    % So I can make y-axis limits all the same for a particular
                    % feature, I need to find out the maximum and minimum
                    % feature values across all participants
                    max_perCond = zeros(numConds, 1);
                    min_perCond = zeros(numConds, 1);
                    for cond = 1:numConds
                        max_perCond(cond) = ...
                            max(plotReady_mean.(condNames{cond})(:, feat_idx, :) + ...
                            plotReady_SEM.(condNames{cond})(:, feat_idx, :), [], 'all');
                        min_perCond(cond) = ...
                            min(plotReady_mean.(condNames{cond})(:, feat_idx, :) - ...
                            plotReady_SEM.(condNames{cond})(:, feat_idx, :), [], 'all');
                    end
                    
                    % Figure out ylim min and max
                    max_yaxis = max(max_perCond);
                    min_yaxis = min(min_perCond);
                    
                    % Loop through conditions (horizontally through subplots)
                    for cond_idx = 1:numConds
                        % Plot in appropriate subplot
                        axes(subs(row_idx*numConds - numConds + cond_idx));
                        
                        % Iterate through condition iterations (i.e., for
                        % neutral, iterate through neutral 1, neutral 2,
                        % ...)
                        hold on
                        for cond_count = 1:conds_perType(cond_idx)
                            % Since conditions have different time lengths,
                            % just create time vectors separately
                            timeVec = 0:(1/Fs):...
                                (size(plotReady_mean.(condNames{cond_idx}), 1)...
                                - 1)/Fs;
                           shadedBounds_linePlot(...
                            plotReady_mean.(condNames{cond_idx})(:, ...
                            feat_idx, cond_count), ...
                            plotReady_mean.(condNames{cond_idx})(:, ...
                            feat_idx, cond_count) + ...
                            plotReady_SEM.(condNames{cond_idx})(:, ...
                            feat_idx, cond_count), ...
                            plotReady_mean.(condNames{cond_idx})(:, ...
                            feat_idx, cond_count) - ...
                            plotReady_SEM.(condNames{cond_idx})(:, ...
                            feat_idx, cond_count), ...
                            timeVec, colorLabels{cond_idx}(cond_count, :), ...
                            markerLabels{cond_count});
                        end
                        
                        % Set ylim
                        ylim([min_yaxis, max_yaxis])
                        
                        % I only want a ylabel on the leftmost plots
                        if cond_idx == 1
                            ylabel(['\Delta ', featNames{feat_idx}]);
                        end
                        % I only want xlabels on the bottom plots
                        if row_idx == featPer || feat_idx == numFeats
                            xlabel('Time (s)');
                        end
                        % I only want subplot titles on the top plots
                        if row_idx == 1
                            title(condNames{cond_idx});
                        end
                        
                        % I want to create one legend per condition type
                        iterNames = cell(conds_perType(cond_idx), 1);
                        for iter = 1:conds_perType(cond_idx)
                            iterNames{iter} = ['Rep. ', num2str(iter)];
                        end
                        if row_idx == 1
                            legend(iterNames);
                        end
                        
                        hold off
                    end
                    
                    % Increment feature index
                    feat_idx = feat_idx + 1;
                end
            end
            
            % Add title to entire subplot grid
            sgtitle([groupName, ' Mean \pm SEM Figure ', num2str(fig_idx), ...
                ' of ', num2str(numFigs)]);
            
            % Save plot
            saveas(fig, [saveDir, groupName, '_meanSEM_', ...
                num2str(fig_idx)], 'fig');
            close all
        end
        
    else
        % Create color map based on number of subjects
        rng(0);     % Set random seed so that colors are always the same
        colorLabels = rand(numSubjects, 3);
        
        % Create cell array of subject IDs based on subject numbers
        subjectIDs = cell(length(subjectNums), 1);
        for sub = 1:numel(subjectIDs)
            subjectIDs{sub} = num2str(subjectNums(sub));
        end
        
        % Start creating figures
        for feat_idx = 1:numFeats
            % Initialize figure
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Figure out how many separate protocol conditions there are
            % (arbitrarily use the first subject)
            numProts = 0;
            for cond_idx = 1:numConds
                numProts = numProts + ...
                    size(featSegm_group{1}.(condNames{cond_idx}), 3);
            end
            
            % Initialize a tight subplot grid
            numCols = 2;
            numRows = ceil(numProts/numCols);
            subs = tight_subplot(numRows, numCols);
            
            % So I can make y-axis limits all the same for a particular
            % feature, I need to find out the maximum and minimum
            % feature values across all participants
            max_perSubj_perCond = zeros(numSubjects, numConds);
            min_perSubj_perCond = zeros(numSubjects, numConds);
            max_perSubj = zeros(numSubjects, 1);
            min_perSubj = zeros(numSubjects, 1);
            for sub = 1:numSubjects
                for cond = 1:numConds
                    max_perSubj_perCond(sub, cond) = max(featSegm_group{...
                        sub}.(condNames{cond})(:, feat_idx, :), [], 'all');
                    min_perSubj_perCond(sub, cond) = min(featSegm_group{...
                        sub}.(condNames{cond})(:, feat_idx, :), [], 'all');
                end
                max_perSubj = max(max_perSubj_perCond, [], 2);
                min_perSubj = min(min_perSubj_perCond, [], 2);
            end
            
            % Figure out ylim min and max
            max_yaxis = max(max_perSubj);
            min_yaxis = min(min_perSubj);
            
            % Loop through rows of subplots in this figure, keeping track
            % of which condition type and exact protocol condition we're on
            cond_idx = 1; cond_count = 1;
            for row_idx = 1:numRows
                for col_idx = 1:numCols
                    % Make sure we still have data to be plotted
                    % (if statement only relevant if number of protocol
                    % conditions, i.e., numProts is odd)
                    if cond_idx <= numConds
                        % Plot in appropriate subplot
                        axes(subs(numCols*(row_idx-1) + col_idx));
                        
                        % Iterate through participants and overlay
                        hold on
                        for sub = 1:numSubjects
                            % The way the data are stored, if a participant does
                            % not have data for a certain protocol condition, then
                            % NaNs are stored. So, we just need to loop through the
                            % subplots and if data are not available, just skip
                            if ~isnan(sum(featSegm_group{sub}.(...
                                    condNames{cond_idx})(:, feat_idx, cond_count)))
                                % Since conditions have different time lengths,
                                % just create time vectors separately
                                timeVec = 0:(1/Fs):...
                                    (size(featSegm_group{sub}.(condNames{cond_idx}), 1)...
                                    - 1)/Fs;
                                plot(timeVec, ...
                                    featSegm_group{sub}.(...
                                    condNames{cond_idx})(:, feat_idx, cond_count), ...
                                    'Color', colorLabels(sub, :));
                            end
                        end
                        hold off
                        
                        % Set ylim
                        ylim([min_yaxis, max_yaxis])
                        
                        % Add subplot title to know which condition this is
                        title([condNames{cond_idx}, ' ', num2str(cond_count)]);
                        
                        % I only want xlabels on the bottom plots
                        if row_idx == numRows || cond_count == ...
                                size(featSegm_group{1}.(condNames{cond_idx}), 3)
                            xlabel('Time (s)');
                        end
                        
                        % I want to create one legend for all, so I will just
                        % arbitrarily assign one of the subplots on the right
                        % to have a legend on the outside bottom right
                        if row_idx == 1 && col_idx == numCols
                            legend(subjectIDs, 'Location', 'southeastoutside')
                        end
                        % I only want a ylabel on the leftmost plots
                        if col_idx == 1
                            ylabel(['\Delta ', featNames{feat_idx}]);
                        end
                        
                        % If we need to reset cond_count and increment cond_idx
                        % (arbitarily use the first subject)
                        if cond_count == ...
                                size(featSegm_group{1}.(condNames{cond_idx}), 3)
                            cond_count = 1;
                            cond_idx = cond_idx + 1;
                        else
                            cond_count = cond_count + 1;
                        end
                        
                    end
                    
                end
            end
            
            % Add title to entire subplot grid
            sgtitle([groupName, ' Subject Overlay Figure for Feature: ', ...
                featNames{feat_idx}]);
            
            % Save plot
            saveas(fig, [saveDir, groupName, '_subjOverlay_', ...
                featNames{feat_idx}], 'fig');
            
            close all
        end
        
    end
    
end

end

