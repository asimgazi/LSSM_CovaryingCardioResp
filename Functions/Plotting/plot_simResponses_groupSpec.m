function plot_simResponses_groupSpec(featPer, titlePrefix, fileNamePrefix, ...
    responses, timeVec, varargin)
% This function takes simulated responses and plots and saves them
% according to inputs provided

% --------- Inputs --------------- %
% featPer: number of feature subplots per figure
% titlePrefix: figure title prefix to use; the suffix will be _1, _2, and 
%              so on based on (# features/featPer)
% fileNamePrefix: filename prefix to use; the suffix will be _1, _2, and so
%                 on based on (# features/featPer)
% responses: cell array of TxMxR arrays, where T is number of time points, M is
%            dimensionality of features/states/etc., and R is the number of
%            simulated responses we have per feature set (i.e., if we
%            simulated responses for 3 different conditions, VNS, trauma,
%            and neutral, then R = 3). Each cell holds one participant's
%            responses
% timeVec: Tx1 time vector to plot against (in seconds)

% (optional)
% 'vtn' or 'vt' or 'vn' or 'v' - If you input one of these strings, you
% are specifying that you want a legend on each plot based on the overlaid
% responses and that the responses will be named according to input cond.
% (e.g., 'vt' - VNS and trauma)

% 'names' - if this flag is provided, the next varargin cell needs to be a
% cell array of names corresponding to the simulated responses provided. If
% no names are provided, the responses are simply labeled by index

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
    if strcmp(varargin{arg}, 'vtn')
        condNames = {'VNS'; 'Trauma'; 'Neutral'};
    elseif strcmp(varargin{arg}, 'vt')
        condNames = {'VNS'; 'Trauma'};
    elseif strcmp(varargin{arg}, 'vn')
        condNames = {'VNS'; 'Neutral'};
    elseif strcmp(varargin{arg}, 'v')
        condNames = {'VNS'};
    elseif strcmp(varargin{arg}, 'names')
        responseNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanSEM = false;    % We will plot responses themselves
        subjectNums = varargin{arg + 1};    % Subject numbers for legend
    end
end

% Set defaults
if ~exist('condNames', 'var')
    % If they don't provide condition names, I'll just number them
    condNames = {};
    for r = 1:size(responses{1}, 3)
        condNames = [condNames; {['Cond. ', num2str(r)]}];
    end
end
if ~exist('responseNames', 'var')
    % If they don't provide response names, I'll just number them
    responseNames = {};
    for m = 1:size(responses{1}, 2)
        responseNames = [responseNames; {['Dim. ', num2str(m)]}];
    end
end
if ~exist('meanSEM', 'var')
    meanSEM = true; % read documentation above for what this bool dictates
end

% Just so we have them as easy to understand variables
numSubjects = numel(responses);
numFeats = size(responses{1}, 2);
numConds = size(responses{1}, 3);

% Find out how many features we're dealing with and figures we'll need
numFigs = ceil(numFeats/featPer);

% Branch off based on meanSEM boolean variable
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
    for cond_idx = 1:numConds
        for sub = 1:numSubjects
            agg_responses{cond_idx} = ...
                cat(3, agg_responses{cond_idx}, ...
                responses{sub}(:, :, cond_idx));
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
    feat_idx = 1;       % To help us iterate through features
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
                    shadedBounds_linePlot(...
                        mean_responses{cond_idx}(:, feat_idx), ...
                        mean_responses{cond_idx}(:, feat_idx) + ...
                        SEM_responses{cond_idx}(:, feat_idx), ...
                        mean_responses{cond_idx}(:, feat_idx) - ...
                        SEM_responses{cond_idx}(:, feat_idx), ...
                        timeVec, colorLabels(cond_idx, :));
                    ylabel(responseNames{feat_idx});
                    xlabel('Time (s)');
                end
                
                % Add legend
                if numel(condNames) > 1
                    legend(condNames);
                end
                hold off
            end
            
            % Increment feature index
            feat_idx = feat_idx + 1;
        end
        
        % Add title to entire subplot grid
        sgtitle([titlePrefix, '- Feat. Set ', num2str(fig_idx)]);
        
        % Save plot
        saveas(fig, [fileNamePrefix, '_', num2str(fig_idx)], ...
            'fig');
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
                max_perSubj = zeros(numSubjects, 1); 
                min_perSubj = zeros(numSubjects, 1);
                for sub = 1:numSubjects
                    max_perSubj(sub) = ...
                        max(responses{sub}(:, feat_idx, :), [], 'all');
                    min_perSubj(sub) = ...
                        min(responses{sub}(:, feat_idx, :), [], 'all');
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
                        plot(timeVec, ...
                            responses{sub}(:, feat_idx, cond_idx), ...
                            'Color', colorLabels(sub, :));
                    end
                    
                    % Set ylim
                    ylim([min_yaxis, max_yaxis])
                    
                    % I only want a ylabel on the leftmost plots
                    if cond_idx == 1
                        ylabel(responseNames{feat_idx});
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
        sgtitle([titlePrefix, '- Feat. Set ', num2str(fig_idx)]);
        
        % Save plot
        saveas(fig, [fileNamePrefix, '_', num2str(fig_idx)], ...
            'fig');
        close all
    end
end


end

