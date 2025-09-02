function [binarySign_forConds, rankOrders_forConds, values_forGroups] = ...
    statTest_LowD_simResponses_groupCompare(dim2plot, plotTitle, fileName, ...
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

% 'groupStat' - indicates you want to statistically test for group
% differences (for closed-loop analyses, i.e., Fig. 7 in TBME paper), 
% rather than for condition differences (Fig. 6 in TBME paper)
% (default is false so default is doing the condition statistics)

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
    elseif strcmp(varargin{arg}, 'groupStat')
        groupStat = true;
        tSteps_ignore = (varargin{arg + 1})*Fs;
        % We need to ignore timesteps for when stimulation is off
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
if ~exist('groupStat', 'var')
    groupStat = false;
end

% Set default outputs
binarySign_forConds = -1;
rankOrders_forConds = -1;
values_forGroups = -1;


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
    

    % Now I can make each column unit variance for PCA transformation
    pcaReady = zeros(size(stackedResps));
    for col = 1:size(stackedResps, 2)
        pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
    end

    % Apply PCA and store away variance explained
    % (outputs I ignore are the the variances
    % in the latent space, Hotelling's T-squared statistic for each
    % observation, and the mean of each column)
    [pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
        'Centered',false);



    % At this point for the statistical analyses, I will need to use the
    % PCA loadings  learned from the mean responses to map each of the
    % participant-specific responses for analyses


    % I am interested in changes along PC1, so I just need the first column
    % of the loadings
    PC1loads = pcaLoads(:, 1);

    % Initialize result
    PC1_groupSep_condSep = cell(numGroups, numConds); % Initialize
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            % First dimension will be time steps for condition
            % Second dimension is number of subjects in group
            PC1_groupSep_condSep{grp_idx, cond_idx} = ...
                zeros(numSamples(cond_idx), numSubjects(grp_idx));
        end
    end

    % Map each of the participant's responses using this PC1 loading
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            for sub = 1:numSubjects(grp_idx)
                PC1_groupSep_condSep{grp_idx, cond_idx}(:, sub) = squeeze(...
                    resps_groupSep_condSep{grp_idx, cond_idx}(:, :, sub)*...
                    PC1loads);
            end
        end
    end

    % Now, we want to compute the mean across time for each condition,
    % still keeping groups separate for now
    PC1means_groupSep = cell(numGroups, 1);
    for grp_idx = 1:numGroups
        % First dimension will be condition ('v', 't', 'n')
        % Second dimension of arrays will remain subject
        PC1means_groupSep{grp_idx} = zeros(numConds, numSubjects(grp_idx));
    end

    % Compute means
    for grp_idx = 1:numGroups
        for cond_idx = 1:numConds
            PC1means_groupSep{grp_idx}(cond_idx, :) = ...
                mean(PC1_groupSep_condSep{grp_idx, cond_idx}, 1);
        end
    end

    % Now we branch depending on if we are comparing groups or conditions
    if groupStat
        % At this point, we just need to compute differences between the
        % just-in-time stimulation and the trauma recall responses for all
        % groups and then return the mean differences in the velocities
        % along PC1

        % Let's first compute velocities
        velocsPC1_groupSep_condSep = cell(numGroups, numConds);
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                % diff() operates along axis 1 by default
                velocsPC1_groupSep_condSep{grp_idx, cond_idx} = ...
                    diff(PC1_groupSep_condSep{grp_idx, cond_idx})./...
                    diff(timeVec(1:size(PC1_groupSep_condSep{grp_idx, cond_idx}, 1)))';
            end
        end

        % Now, I need to compute the difference between velocities for the
        % first two conditions for each group separately (I know
        % just-in-time stimulation is condition 1 and trauma recall is
        % condition 2)
        velocDiffs_PC1_groupSep = cell(numGroups, 1);
        for grp_idx = 1:numGroups
            velocDiffs_PC1_groupSep{grp_idx} = ...
                velocsPC1_groupSep_condSep{grp_idx, 1} - ...
                velocsPC1_groupSep_condSep{grp_idx, 2};
        end

        % Now I compute the mean velocity difference for each subject
        % But we only care about differences once stimulation is turned on
        % Note that with velocities, we have one less timestep, so I just
        % start from tSteps_ignore rather than tSteps_ignore + 1
        mean_velocDiffs_PC1_groupSep = cell(numGroups, 1);
        for grp_idx = 1:numGroups
            mean_velocDiffs_PC1_groupSep{grp_idx} = mean(...
                velocDiffs_PC1_groupSep{grp_idx}(tSteps_ignore:end, :), 1);
        end

        % Assign and ready to return
        % What is returned is a cell array with length of number of groups
        % Each cell contains a 1-D array with a mean velocity difference
        % between condition 1 and 2 for each subject
        values_forGroups = mean_velocDiffs_PC1_groupSep;


    else
        % We do not care about groups in this case and just want to show
        % that the condition responses go in the direction wanted

        % I will need to output a 1 or a 0 for each participant for each
        % condition indicating whether the mean response was positive or
        % negative (1 for positive and 0 for negative)

        % First dimension will be subject and second dimension will be
        % condition
        binarySign_forConds = -1*ones(sum(numSubjects), numConds);

        % Check sign and assign 0 or 1 accordingly
        for grp_idx = 1:numGroups
            for cond_idx = 1:numConds
                for sub = 1:numSubjects(grp_idx)
                    if PC1means_groupSep{grp_idx}(cond_idx, sub) < 0
                        % To iterate properly, we need to factor in
                        % subjects from previous groups if needed
                        if grp_idx > 1
                            % Value is negative so assign 0
                            binarySign_forConds(sum(numSubjects(1:grp_idx-1)) + sub, cond_idx) = 0;
                        else
                            binarySign_forConds(sub, cond_idx) = 0;
                        end
                    else
                        % To iterate properly, we need to factor in
                        % subjects from previous groups if needed
                        if grp_idx > 1
                            % Value is positive so assign 1
                            binarySign_forConds(sum(numSubjects(1:grp_idx-1)) + sub, cond_idx) = 1;
                        else
                            binarySign_forConds(sub, cond_idx) = 1;
                        end
                    end
                end
            end
        end

        %  I will also need to output the rank ordering per participant of
        % stimulation, trauma, and neutral so that I can then do a repeated
        % measures correlation between simulated and real

        % First dimension will be subjects and second dimension will be
        % condition
        rankOrders_forConds = -1*ones(sum(numSubjects), numConds);

        % What I will do is rank the values from 1 to 3, where 3 is
        % greatest value across conditions, 2 is middle, and 1 is smallest
        for grp_idx = 1:numGroups
            for sub = 1:numSubjects(grp_idx)
                % Go through all six possible ranking combinations
                if PC1means_groupSep{grp_idx}(1, sub) < ...
                        PC1means_groupSep{grp_idx}(2, sub) && ...
                        PC1means_groupSep{grp_idx}(2, sub) < ...
                        PC1means_groupSep{grp_idx}(3, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [1, 2, 3];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [1, 2, 3];
                    end

                elseif PC1means_groupSep{grp_idx}(1, sub) > ...
                        PC1means_groupSep{grp_idx}(2, sub) && ...
                        PC1means_groupSep{grp_idx}(2, sub) > ...
                        PC1means_groupSep{grp_idx}(3, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [3, 2, 1];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [3, 2, 1];
                    end

                elseif PC1means_groupSep{grp_idx}(2, sub) < ...
                        PC1means_groupSep{grp_idx}(1, sub) && ...
                        PC1means_groupSep{grp_idx}(1, sub) < ...
                        PC1means_groupSep{grp_idx}(3, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [2, 1, 3];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [2, 1, 3];
                    end

                elseif PC1means_groupSep{grp_idx}(2, sub) > ...
                        PC1means_groupSep{grp_idx}(1, sub) && ...
                        PC1means_groupSep{grp_idx}(1, sub) > ...
                        PC1means_groupSep{grp_idx}(3, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [2, 3, 1];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [2, 3, 1];
                    end

                elseif PC1means_groupSep{grp_idx}(1, sub) < ...
                        PC1means_groupSep{grp_idx}(3, sub) && ...
                        PC1means_groupSep{grp_idx}(3, sub) < ...
                        PC1means_groupSep{grp_idx}(2, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [1, 3, 2];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [1, 3, 2];
                    end

                elseif PC1means_groupSep{grp_idx}(2, sub) < ...
                        PC1means_groupSep{grp_idx}(3, sub) && ...
                        PC1means_groupSep{grp_idx}(3, sub) < ...
                        PC1means_groupSep{grp_idx}(1, sub)
                    if grp_idx > 1
                        % Value is negative so assign 0
                        rankOrders_forConds(sum(numSubjects(1:grp_idx-1)) + sub, :) = ...
                            [3, 1, 2];
                    else
                        rankOrders_forConds(sub, :) = ...
                            [3, 1, 2];
                    end
                end
                
            end
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
        video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
            [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
            'mrkrNames', groupNames, 'clrNames', condNames, ...
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
        video_lowDtrajectories(plotReady, timeVec, dim2plot, ...
            [fileName, ' - Overlaid'], 'plotTitle', plotTitle, ...
            'mrkrNames', groupNames, 'clrNames', condNames, ...
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

