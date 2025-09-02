function plotLowD_simResponses_groupCompare(dim2plot, plotTitle, fileName, ...
    responses_groupSep, Fs, varargin)
% This function will take simulated responses of different groups of data,
% stack them all up, apply PCA, and then reseparate according to groups
% according to inputs provided

% ---------- Inputs -------------------- %
% dim2plot: scalar = 2 or 3 that dictates whether to plot in 2-D or 3-D
% plotTitle: figure title to use
% fileName: filename to use
% responses_groupSep: cell array of cell arrays. First layer of cells
%                     separates groups. Second layer of cells separates
%                     participants within the groups. Each cell contains a 
%                     struct with R fields. Each field contains a T_i x M
%                     array. T_i represents the number of timepoints for
%                     the ith response, i \in {1, ..., R}. M is the number
%                     of features (number of features is same across fields
%                     (e.g., if we simulated responses for 3 different 
%                     conditions, VNS, trauma, and neutral, then R = 3)
%                     (e.g., if we have 60 samples of trauma data, 60
%                     samples of neutral data, and 120 samples of VNS data,
%                     then T_trauma = 60, T_neutral = 60, and T_VNS = 120)
%                     NOTE: the struct fieldnames are used as condition names
%                     unless they are input separately (see 'condNames'
%                     flag)
% Fs: sampling rate in Hz

% (optional)
% 'verbose' - if you want to see all warnings (default is not verbose)

% 'subjOverlay' - if this flag is provided, the next varargin cell needs to
% be a cell array with a vector of subject numbers for each group in each 
% cell to use for legends. What this flag entails is it will set 
% meanResp = false. This will overlay all subjects' responses from a group for a 
% particular condition on top of each other while separating groups and
% conditions into separate plots. If meanResp = true, means will be computed across 
% participants for each group and condition and then overlaid on the same plot for 
% easy comparison, so the legend will just be based on condition and group

% 'groupNames'- will need to be followed by a cell array of groups names
% corresponding to responses_groupSep

% 'condNames' - will need to be followed by a cell array of condition names
% corresponding to responses_groupSep (dimension R)

% 'additional' - flag to inform the function that additional responses will
% be included beyond what is included in responses_groupSep. These
% additional responses will follow the same form as
% responses_groupSep (M needs to match; T_i and R can differ). The
% additional responses need to be input after this flag. These additional
% responses will not be used when learning the PCA transformation matrix,
% but will be transformed into same subspace as other responses and plotted

% 'add_groupNames' - serves same purpose as 'groupNames' flag but for the
% additional responses

% 'add_condNames' - serves same purpose as 'condNames' flag but for the
% additional responses

% 'scaleSeparate' - flag that specifies boolean that responses need to be
% scaled for unit variance separately before being linearly transformed
% using PCA transformation matrix; otherwise, the variances of the first
% set of responses are used to scale the additional responses too

% concatAdditional' - flag that specifies boolean that additional responses
% should be concatenated with primary responses when applying PCA such that
% both sets of data are considered when scaling to unit variance and
% learning PCA transformation matrix; otherwise, additional responses are
% not considered when applying PCA

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'condNames')
        condNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'verbose')
        verbose = true;
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanResp = false;    % We will plot responses themselves
        subjectNums = varargin{arg + 1};    % Subject numbers for legend
    elseif strcmp(varargin{arg}, 'groupNames')
        groupNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'additional')
        addResps = true;    % Let the function know there are additionals
        addResps_groupSep = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'add_groupNames')
        add_groupNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'add_condNames')
        add_condNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'scaleSeparate')
        scaleSeparate = true;
    elseif strcmp(varargin{arg}, 'concatAdditional')
        concatAdditional = true;
    end
end

% Set defaults
if ~exist('condNames', 'var')
    % If they don't provide condition names use the input structs' field
    % names
    % (arbitrarily used first group's first subject)
    condNames = fieldnames(responses_groupSep{1}{1});
end
if ~exist('groupNames', 'var')
    % If they don't provide condition names, I'll just number them
    groupNames = {};
    for g = 1:numel(responses_groupSep)
        groupNames = [groupNames; {['Group ', num2str(g)]}]; 
    end
end
if ~exist('verbose', 'var')
    verbose = false;
end
if ~exist('meanResp', 'var')
    meanResp = true; % read documentation above for what this bool dictates
end

% If we have additional responses...
if ~exist('addResps', 'var')
    addResps = false;
end
if addResps
    if ~exist('add_condNames', 'var')
        % If they don't provide condition names for the additional
        % responses, use the structs' field names
        % (arbitrarily used first group's first subject)
        add_condNames = fieldnames(addResps_groupSep{1}{1});
    end
    if ~exist('add_groupNames', 'var')
        % If they don't provide condition names, I'll just number them
        add_groupNames = {};
        for g = 1:numel(addResps_groupSep)
            add_groupNames = [add_groupNames; {['Add. Group ', num2str(g)]}];
        end
    end
    if ~exist('scaleSeparate', 'var')
        scaleSeparate = false;
    end
    if ~exist('concatAdditional', 'var')
        concatAdditional = false;
    end
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
numGroups = numel(responses_groupSep);
numSubjects = zeros(numGroups, 1);
for grp_idx = 1:numGroups
    numSubjects(grp_idx) = numel(responses_groupSep{grp_idx});
end
numConds = length(condNames);
% Arbitrarily use the first group's first subject
Rnames = fieldnames(responses_groupSep{1}{1});  % Separate in case condNames differ
numSamples = zeros(numConds, 1);
for cond_idx = 1:numConds
    % Arbitrarily use the first group's first subject
    numSamples(cond_idx) = size(responses_groupSep{1}{1}.(Rnames{cond_idx}), 1);
end

% Arbitrarily use first group's first subject's first field's data
numFeats = size(responses_groupSep{1}{1}.(Rnames{1}), 2);

% Save these away if additional responses are provided
if addResps
    add_numGroups = numel(addResps_groupSep);
    add_numSubjects = zeros(add_numGroups, 1);
    for grp_idx = 1:add_numGroups
        add_numSubjects(grp_idx) = numel(addResps_groupSep{grp_idx});
    end
    add_numConds = length(add_condNames);
    add_Rnames = fieldnames(addResps_groupSep{1}{1}); % Separate in case add_condNames differ
    add_numSamples = zeros(add_numConds, 1);
    for cond_idx = 1:add_numConds
        % Arbitrarily use the first group's first subject
        add_numSamples(cond_idx) = size(addResps_groupSep{1}{1}.(add_Rnames{cond_idx}), 1);
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
                responses_groupSep{grp_idx}{sub}.(Rnames{cond_idx});
        end
    end
end


% Do the same separation and filling for the additional responses if they
% are provided
if addResps
    addResps_groupSep_condSep = cell(add_numGroups, add_numConds);
    for grp_idx = 1:add_numGroups
        for cond_idx = 1:add_numConds
            % Initialize
            addResps_groupSep_condSep{grp_idx, cond_idx} = ...
                zeros(add_numSamples(cond_idx), numFeats, add_numSubjects(grp_idx));
        end
    end
    % Now fill the cell array's arrays
    for grp_idx = 1:add_numGroups
        for cond_idx = 1:add_numConds
            for sub = 1:add_numSubjects(grp_idx)
                addResps_groupSep_condSep{grp_idx, cond_idx}(:, :, sub) = ...
                    addResps_groupSep{grp_idx}{sub}.(add_Rnames{cond_idx});
            end
        end
    end
end


% To encode time, I am going to create a time vector according to the
% longest time protocol condition we have and then just use a subset of
% that for the other conditions
if addResps
    timeVec = 0:(1/Fs):(max([numSamples; add_numSamples])-1)/Fs;
else
    timeVec = 0:(1/Fs):(max(numSamples)-1)/Fs;
end

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
    grpLabels = {'o', 'p', 's', '^', 'd', '*', '.', 'x', '_', '|', ...
        '+', 'v', '>', '<', 'h'};
    
    % Compute means within each condition, for each group
    % (3rd dimension is subject)
    meanResps_groupSep_condAvg = cell(numGroups, numConds);
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            meanResps_groupSep_condAvg{grp_idx, cond_idx} = ...
                mean(resps_groupSep_condSep{grp_idx, cond_idx}, 3);
        end
    end
    % Do the same averaging for the additional responses if provided
    if addResps
        add_meanResps_groupSep_condAvg = cell(add_numGroups, add_numConds);
        for grp_idx = 1:add_numGroups
            for cond_idx = 1:add_numConds
               add_meanResps_groupSep_condAvg{grp_idx, cond_idx} = ...
                    mean(addResps_groupSep_condSep{grp_idx, cond_idx}, 3);
            end
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
                meanResps_groupSep_condAvg{grp_idx, cond_idx}];
            groupLabels = [groupLabels; grp_idx*...
                ones(size(meanResps_groupSep_condAvg{grp_idx, cond_idx}, 1), 1)];
            condLabels = [condLabels; cond_idx*...
                ones(size(meanResps_groupSep_condAvg{grp_idx, cond_idx}, 1), 1)];
        end
    end
    
    % Check to make sure sizes are correct
    if size(stackedResps, 1) ~= sum(numSamples)*numGroups || ...
            length(condLabels) ~= sum(numSamples)*numGroups || ...
            length(groupLabels) ~= sum(numSamples)*numGroups
        error('Stacking did not go as planned');
    end
    % Perform the same stacking for the additional responses if provided
    if addResps
        % Initialize stacked responses
        add_stackedResps = [];
        add_groupLabels = [];
        add_condLabels = [];
        
        % Loop and stack
        for grp_idx = 1:add_numGroups
            for cond_idx = 1:add_numConds
                add_stackedResps = [add_stackedResps; ...
                    add_meanResps_groupSep_condAvg{grp_idx, cond_idx}];
                add_groupLabels = [add_groupLabels; grp_idx*...
                    ones(size(add_meanResps_groupSep_condAvg{grp_idx, cond_idx}, 1), 1)];
                add_condLabels = [add_condLabels; cond_idx*...
                    ones(size(add_meanResps_groupSep_condAvg{grp_idx, cond_idx}, 1), 1)];
            end
        end
        
        % Check to make sure sizes are correct
        if size(add_stackedResps, 1) ~= sum(add_numSamples)*add_numGroups || ...
                length(add_condLabels) ~= sum(add_numSamples)*add_numGroups || ...
                length(add_groupLabels) ~= sum(add_numSamples)*add_numGroups
            error('Stacking did not go as planned');
        end
    end
    
    
    % If we need to combine the primary responses with the additional
    % responses for PCA purposes, do that now and then reseparate
    if addResps && concatAdditional
        % Keep track of how many rows are from primary responses
        numRows_primary = size(stackedResps, 1);
        
        % Stack them on top of each other
        comboStack = [stackedResps; add_stackedResps];
        
        % Make each column unit variance
        pcaReady = zeros(size(comboStack));
        for col = 1:size(comboStack, 2)
            pcaReady(:, col) = (1/std(comboStack(:, col)))*comboStack(:, col);
        end
        
        % Apply PCA and store away variance explained and the scores (i.e.,
        % input matrix transformed into principal component space)
        [pcaLoads, combo_pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % Separate back into primary responses and additional responses
        pcaResult = combo_pcaResult(1:numRows_primary, :);
        add_pcaResult = combo_pcaResult((numRows_primary + 1):end, :);
    else
        % Otherwise, apply PCA only to the primary responses and then
        % transform the additional responses separately
        
        % Now I can make each column unit variance for PCA transformation
        pcaReady = zeros(size(stackedResps));
        for col = 1:size(stackedResps, 2)
            pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
        end
        
        % If additional responses are provided, scale according to user prefer
        if addResps
            add_pcaReady = zeros(size(add_stackedResps));
            if scaleSeparate
                % If we scale separately, we use add_stackedResps SDs to scale
                for col = 1:size(add_stackedResps, 2)
                    add_pcaReady(:, col) = (1/std(add_stackedResps(:, col)))*add_stackedResps(:, col);
                end
            else
                % If we scale same, we use stackedResps SDs to scale
                for col = 1:size(add_stackedResps, 2)
                    add_pcaReady(:, col) = (1/std(stackedResps(:, col)))*add_stackedResps(:, col);
                end
            end
        end
        
        % Apply PCA and store away variance explained
        % (outputs I ignore are the the variances
        % in the latent space, Hotelling's T-squared statistic for each
        % observation, and the mean of each column)
        [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);
        
        % Transform additional responses using the PCA transformation matrix
        if addResps
            add_pcaResult = add_pcaReady*pcaLoads;
        end
    end
    
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
    
    % Separate the additional responses as well if provided
    if addResps
        add_plotReady = cell(add_numGroups, add_numConds);
        for grp_idx = 1:add_numGroups
            for cond_idx = 1:add_numConds
                add_plotReady{grp_idx, cond_idx} = ...
                    add_pcaResult(add_groupLabels == grp_idx & ...
                    add_condLabels == cond_idx, :);
                
                % Check size to make sure everything was stored properly
                if size(add_plotReady{grp_idx, cond_idx}, 1) ~= add_numSamples(cond_idx) || ...
                        size(add_plotReady{grp_idx, cond_idx}, 2) ~= numFeats
                    error('add_plotReady sizing is incorrect')
                end
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
            if addResps && concatAdditional
                xlabel(['PC1 primary and add. data var. expl. : ', ...
                    num2str(varExplain(1)), '%']);
                ylabel(['PC2 primary and add. data var. expl. : ', ...
                    num2str(varExplain(2)), '%']);
            else
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
            end
            
            % Include legend
            if numel(condNames) > 1
                legend(condNames, 'Location', 'best')
                legend('boxoff')
            end
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
        if addResps && concatAdditional
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames, ...
                'xlabel', ...
                ['PC1 primary and add. data var. expl. : ', ...
                num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 primary and add. data var. expl. : ', ...
                num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 primary and add. data var. expl. : ', ...
                num2str(varExplain(3)), '%'], 'velocities');
        else
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
        end
        
        
        % Close all figures
        close all
        
        
        % Create the same plots for the additional responses if provided
        if addResps
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, add_numGroups);
            
            % Iterate through groups and conditions and overlay
            for grp_idx = 1:add_numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:add_numConds
                    % Plot one scatter point at a time to leverage time labels
                    for t = 1:size(add_plotReady{grp_idx, cond_idx}, 1)
                        s = scatter(add_plotReady{grp_idx, cond_idx}(t, 1), ...
                            add_plotReady{grp_idx, cond_idx}(t, 2), sizeLabels(t), ...
                            colorLabels(cond_idx, :), grpLabels{grp_idx}, 'filled');
                        s.MarkerFaceAlpha = alphaLabels(t);
                        
                        % If this is the final datapoint, then store legend
                        % info
                        if t == size(add_plotReady{grp_idx, cond_idx}, 1)
                            s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                        else
                            s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                if concatAdditional
                    xlabel(['PC1 primary and add. data var. expl. : ', ...
                        num2str(varExplain(1)), '%']);
                    ylabel(['PC2 primary and add. data var. expl. : ', ...
                        num2str(varExplain(2)), '%']);
                else
                    xlabel('Projection Along PC1 of Other Data');
                    ylabel('Projection Along PC2 of Other Data');
                end
                
                
                % Include legend
                if numel(add_condNames) > 1
                    legend(add_condNames, 'Location', 'best')
                    legend('boxoff')
                end
                hold off
                
                % Add title to subplot
                title(add_groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' Additional Responses - Separate']);
            
            % Save plot
            saveas(fig, [fileName, '_addResps - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Create video for all data overlaid
            if concatAdditional
                video_lowDtrajectories(add_plotReady, timeVec, dim2plot, ...
                    [fileName, '_addResps - Overlaid'], 'plotTitle', ...
                    [plotTitle, ' Additional Responses'], ...
                    'mrkrNames', add_groupNames, 'clrNames', add_condNames, ...
                    'xlabel', ...
                    ['PC1 primary and add. data var. expl. : ', ...
                    num2str(varExplain(1)), '%'], ...
                    'ylabel', ...
                    ['PC2 primary and add. data var. expl. : ', ...
                    num2str(varExplain(2)), '%'], ...
                    'zlabel', ...
                    ['PC3 primary and add. data var. expl. : ', ...
                    num2str(varExplain(3)), '%'], 'velocities');
            else
                video_lowDtrajectories(add_plotReady, timeVec, dim2plot, ...
                    [fileName, '_addResps - Overlaid'], 'plotTitle', ...
                    [plotTitle, ' Additional Responses'], ...
                    'mrkrNames', add_groupNames, 'clrNames', add_condNames, ...
                    'xlabel', ...
                    'Projection Along PC1 of Other Data', ...
                    'ylabel', ...
                    'Projection Along PC2 of Other Data', ...
                    'zlabel', ...
                    'Projection Along PC3 of Other Data', 'velocities');
            end
            
            % Close all figures
            close all
        end
        
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
            if addResps && concatAdditional
                xlabel(['PC1 primary and add. data var. expl. : ', ...
                    num2str(varExplain(1)), '%']);
                ylabel(['PC2 primary and add. data var. expl. : ', ...
                    num2str(varExplain(2)), '%']);
                zlabel(['PC3 primary and add. data var. expl. : ', ...
                    num2str(varExplain(3)), '%']);
            else
                xlabel(['PC1 Var. Expl. : ', num2str(varExplain(1)), '%']);
                ylabel(['PC2 Var. Expl. : ', num2str(varExplain(2)), '%']);
                zlabel(['PC3 Var. Expl. : ', num2str(varExplain(3)), '%']);
            end
            
            % Include legend
            legend(condNames, 'Location', 'best')
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
        if addResps && concatAdditional
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames, ...
                'xlabel', ...
                ['PC1 primary and add. data var. expl. : ', ...
                num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 primary and add. data var. expl. : ', ...
                num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 primary and add. data var. expl. : ', ...
                num2str(varExplain(3)), '%'], 'velocities');
        else
            video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
                [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
                'mrkrNames', groupNames, 'clrNames', condNames, ...
                'xlabel', ...
                ['PC1 Var. Expl. : ', num2str(varExplain(1)), '%'], ...
                'ylabel', ...
                ['PC2 Var. Expl. : ', num2str(varExplain(2)), '%'], ...
                'zlabel', ...
                ['PC3 Var. Expl. : ', num2str(varExplain(3)), '%'], 'velocities');
        end
        
        % Close all figures
        close all
        
        
        % Create the same plots for the additional responses if provided
        if addResps
            % Create figure for separate subplot for each group
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize subplot grid, one subplot per group
            subs = tight_subplot(1, add_numGroups);
            
            % Iterate through groups and conditions and overlay
            for grp_idx = 1:add_numGroups
                axes(subs(grp_idx));
                hold on
                for cond_idx = 1:add_numConds
                    % Plot one scatter point at a time to leverage time labels
                    for t = 1:size(add_plotReady{grp_idx, cond_idx}, 1)
                        s = scatter3(add_plotReady{grp_idx, cond_idx}(t, 1), ...
                            add_plotReady{grp_idx, cond_idx}(t, 2), ...
                            add_plotReady{grp_idx, cond_idx}(t, 3), ...
                            sizeLabels(t), ...
                            colorLabels(cond_idx, :), grpLabels{grp_idx}, 'filled');
                        s.MarkerFaceAlpha = alphaLabels(t);
                        
                        % If this is the final datapoint, then store legend
                        % info
                        if t == size(add_plotReady{grp_idx, cond_idx}, 1)
                            s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                        else
                            s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        end
                    end
                end
                % Label axes and title
                % Include variance explained by each of the principal components in the
                % axis labels
                if concatAdditional
                    xlabel(['PC1 primary and add. data var. expl. : ', ...
                        num2str(varExplain(1)), '%']);
                    ylabel(['PC2 primary and add. data var. expl. : ', ...
                        num2str(varExplain(2)), '%']);
                    zlabel(['PC3 primary and add. data var. expl. : ', ...
                        num2str(varExplain(3)), '%']);
                else
                    xlabel('Projection Along PC1 of Other Data');
                    ylabel('Projection Along PC2 of Other Data');
                    zlabel('Projection Along PC3 of Other Data');
                end
                
                % Include legend
                legend(add_condNames, 'Location', 'best')
                legend('boxoff')
                hold off
                
                % Add title to subplot
                title(add_groupNames{grp_idx});
            end
            
            % Add title to entire figure
            sgtitle([plotTitle, ' Additional Responses - Separate']);
            
            % Save plot
            saveas(fig, [fileName, '_addResps - Separate'], 'fig');
            
            % Close figures
            close all
            
            % Create video for all data overlaid
            if concatAdditional
                video_lowDtrajectories(add_plotReady, timeVec, dim2plot, ...
                    [fileName, '_addResps - Overlaid'], 'plotTitle', ...
                    [plotTitle, ' Additional Responses'], ...
                    'mrkrNames', add_groupNames, 'clrNames', add_condNames, ...
                    'xlabel', ...
                    ['PC1 primary and add. data var. expl. : ', ...
                    num2str(varExplain(1)), '%'], ...
                    'ylabel', ...
                    ['PC2 primary and add. data var. expl. : ', ...
                    num2str(varExplain(2)), '%'], ...
                    'zlabel', ...
                    ['PC3 primary and add. data var. expl. : ', ...
                    num2str(varExplain(3)), '%'], 'velocities');
            else
                video_lowDtrajectories(add_plotReady, timeVec, dim2plot, ...
                    [fileName, '_addResps - Overlaid'], 'plotTitle', ...
                    [plotTitle, ' Additional Responses'], ...
                    'mrkrNames', add_groupNames, 'clrNames', add_condNames, ...
                    'xlabel', ...
                    'Projection Along PC1 of Other Data', ...
                    'ylabel', ...
                    'Projection Along PC2 of Other Data', ...
                    'zlabel', ...
                    'Projection Along PC3 of Other Data', 'velocities');
            end
            
            % Close all figures
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
    
else
    if addResps
        error(['Function not yet ready for additional responses with ', ...
            'overlaid subjects']);
    end
        
    % Create color map based on total number of subjects
    rng(0);     % Set random seed so that colors are always the same
    colorLabels = rand(sum(numSubjects), 3);
    
    % To differentiate groups, we will use different markers
    grpLabels = {'o', 'p', 's', '^', 'd', '*', '.', 'x', '_', '|', ...
        '+', 'v', '>', '<', 'h'};
    
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
    
    % Since I am not doing any averaging, I essentially need to stack and
    % flatten the responses to be a 2-D array of dimension
    % sum(numSubjects)*sum(numSamples) x numFeats. At the 
    % same time, I will need to keep track of groupLabels, subjLabels and 
    % condLabels so I can later separate when plotting
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
                if size(plotReady{grp_idx, cond_idx}{sub_idx}, 1) ~= numSamples(cond_idx) || ...
                        size(plotReady{grp_idx, cond_idx}{sub_idx}, 2) ~= numFeats
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
                            sizeLabels(t), ...
                            colorLabels(sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
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
                title([groupNames{grp_idx}, '-', condNames{cond_idx}]);
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
                            sizeLabels(t), ...
                            colorLabels(sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
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
            title(condNames{cond_idx});
            
            hold off
        end
        
        % Add title to entire figure
        sgtitle([plotTitle, ' - Overlaid']);
        
        % Save plot
        saveas(fig, [fileName, ' - Overlaid'], 'fig');
        
        % Close both figures
        close all
        
    else
        % Create figure
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
                            sizeLabels(t), ...
                            colorLabels(sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
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
                
                % Set xlim, ylim, and zlim
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
                    legend(subjectIDs{grp_idx}, 'Location', 'eastoutside')
                end
                
                % Title it according to group and condition
                title([groupNames{grp_idx}, '-', condNames{cond_idx}]);
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
                            sizeLabels(t), ...
                            colorLabels(sum(numSubjects(1:(grp_idx-1))) + sub_idx, ...
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
            title(condNames{cond_idx});
            
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

end

