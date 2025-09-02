function plot_realResponses_subjSpec(subjectID, saveDir, featSegm, ...
    Fs, varargin)
% This function takes real responses and plots and saves them
% according to inputs provided

% --------- Inputs --------------- %
% subjectID: string specifying subject ID to use when saving
% saveDir: string specifying the directory where plots should be saved
% featSegm: struct that contains a field for each protocol condition, where
%           within each field, an array exists that is TxMxC, where T is
%           the number of time points in that condition segment, M is the
%           number of features, and C is the number of condition
%           repetitions for that condition
% Fs: sampling frequency; used to construct time vectors to plot against
% (optional)
% 'plotOption': flag that specifies to look for the option for plotting in
%               next function input. Options to choose from:
%               'condSep': each protocol condition separately (i.e.,
%               neutral 1 is plotted separate from neutral 2, etc.)
%               'condAvg': if this option is selected, then the next input
%               should be the number of features per plot since conditions
%               are overlaid on same subplot
%               'condAvg_train': if this option is selected, then the next
%               two inputs should be number of features and the train-test
%               split condition code (e.g., 6.5 for split between 6 and 7)
% 'names': flag that specifies next input as cell array of feature names.
%          If no names are provided, the responses are labeled by index


% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'condSep')
        condAvg = false;
        trainOnly = false;
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
    for m = 1:size(featSegm.n, 2)
        featNames = [featNames; {['Dim. ', num2str(m)]}];
    end
end

% Just so we have them as easy to understand variables
numFeats = length(featNames);
numConds = length(fieldnames(featSegm));

% Create basic color labels (blue, red, green) for protocol conditions
colorLabels = (1/255)*[221, 85, 255; ...
    255, 85, 85; ...
    95, 211, 188];

if condAvg
    % Find out how many features we're dealing with and figures we'll need
    numFigs = ceil(numFeats/featPer);
    
    % Initialize a condition-averaged response struct to store result
    feat_condAvg = struct();
    condNames = fieldnames(featSegm);
    for cond = 1:numConds
        % Initialized to NaNs as default in case no data are available to
        % replace (will ignore NaN responses during plotting)
        feat_condAvg.(condNames{cond}) = ...
            nan*ones(size(featSegm.(condNames{cond}), 1), numFeats);
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
        for cond = 1:numConds
            % Average along the 3rd dimension to average over repetitions
            % Only average the training data
            feat_condAvg.(condNames{cond}) = ...
                mean(featSegm.(condNames{cond})(:, :, 1:trainNums(cond)), ...
                3, 'omitnan');
        end
    else
        % In this case, we can just average over whatever we have since no
        % train/test considerations
        for cond = 1:numConds
            % Average along the 3rd dimension to average over repetitions
            feat_condAvg.(condNames{cond}) = ...
                mean(featSegm.(condNames{cond}), 3, 'omitnan');
        end
    end
    
    % Now that we have the condition averages for this subject, we can plot
    feat_idx = 1;           % To help us iterate through features
    for fig_idx = 1:numFigs
        % Initialize
        fig = figure(1);
        set(fig,'Visible','on');
        
        % Initialize a tight subplot
        subs = tight_subplot(featPer, 1);
        
        % Iterate through subplots
        for subP_idx = 1:featPer
            % Check to make sure we still have features left to plot
            if feat_idx <= numFeats
                % Plot in subplot
                axes(subs(subP_idx));
                
                % To store legend
                legendNames = {};
                
                % Iterate through conditions and overlay
                hold on
                for cond = 1:numConds
                    % Make sure this response is not just all NaNs
                    if ~isnan(sum(feat_condAvg.(condNames{cond})(:, ...
                            feat_idx)))
                        % Create a time vector for this condition
                        timeVec = 0:(1/Fs):...
                            (length(feat_condAvg.(condNames{cond})(:, ...
                            feat_idx))-1)/Fs;
                        plot(timeVec, ...
                            feat_condAvg.(condNames{cond})(:, feat_idx), ...
                            'Color', colorLabels(cond, :), 'LineWidth', 2)
                        ylabel(['\Delta ', featNames{feat_idx}]);
                        xlabel('Time (s)');
                        
                        % Add this condition to the legend since we have
                        % data
                        legendNames{cond} = condNames{cond};
                    end
                end
                
                % Add legend
                legend(legendNames);
                hold off
            end
            
            % Increment feature index
            feat_idx = feat_idx + 1;
        end
        
        % Add title to entire subplot grid
        sgtitle([subjectID, ' Figure ', num2str(fig_idx), ...
            ' of ', num2str(numFigs)]);
        
        % Save plot
        saveas(fig, [saveDir, 'Subject ', subjectID, '_', ...
            num2str(fig_idx)], 'fig');
        close all
    end
else
    % If we do not want to average over repetitions of a condition, then we
    % will be creating a plot for each feature in its own folder
    
    % Figure out how many total plots to create by ignoring NaN repetitions
    totalPlots = 0;
    condNames = fieldnames(featSegm);
    for cond = 1:numConds
        numReps = size(featSegm.(condNames{cond}), 3);
        for rep = 1:numReps
            if ~isnan(sum(featSegm.(condNames{cond})(:, :, rep), 'all'))
                totalPlots = totalPlots + 1;
            end
        end
    end
    
    % We will have two columns of subplots, so figure out how many rows so
    % we know how we can initialize the subplot figures
    numRows = ceil(totalPlots/2);
    
    % Loop through features, creating a figure per feature
    for feat_idx = 1:numFeats
        fig = figure(1);
        set(fig,'Visible','on');
        
        % Initialize a tight subplot
        subs = tight_subplot(numRows, 2);
        
        % Loop through conditions and repetitions, keeping track of subplot
        % index
        plot_idx = 1;
        for cond = 1:numConds
            % Number of repetitions
            numReps = size(featSegm.(condNames{cond}), 3);
            for rep = 1:numReps
                % If data are available to plot
                if ~isnan(sum(featSegm.(condNames{cond})(:, :, rep), 'all'))
                    % Plot in subplot
                    axes(subs(plot_idx));
                    
                    % Create a time vector
                    timeVec = 0:(1/Fs):...
                        (length(featSegm.(condNames{cond})(:, feat_idx, rep))-1)/Fs;
                    
                    % Plot the segmented feature vs. time
                    % We color code based on condition
                    plot(timeVec, ...
                        featSegm.(condNames{cond})(:, feat_idx, rep), ...
                        'Color', colorLabels(cond, :), 'LineWidth', 2)
                    legend(condNames{cond});
                    
                    % Only include ylabels for left plots
                    if mod(plot_idx, 2) == 1
                        ylabel(['\Delta ', featNames{feat_idx}]);
                    end
                    
                    % Only include xlabels for bottom plots
                    if numel(subs) - plot_idx < 2
                        xlabel('Time (s)');
                    end
                    
                    % Increment plot index
                    plot_idx = plot_idx + 1;
                end
            end
        end
        
        % Add title to entire subplot grid
        sgtitle([subjectID, ' ', featNames{feat_idx}]);
        
        % Check if directory exists for this feature; if not, create
        if ~exist([saveDir, featNames{feat_idx}], 'dir')
            mkdir([saveDir, featNames{feat_idx}])
        end
        
        % Save plot
        saveas(fig, [saveDir, featNames{feat_idx}, filesep, subjectID], ...
            'fig');
        close all
    end
    
end

end

