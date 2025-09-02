function plotLowD_simResponses_groupSpec(dim2plot, plotTitle, fileName, ...
    responses, timeVec, varargin)
% This function will take simulated responses, apply PCA to them, and then
% plot and save them according to inputs provided

% ---------- Inputs -------------------- %
% dim2plot: scalar = 2 or 3 that dictates whether to plot in 2-D or 3-D
% plotTitle: figure title to use
% fileName: filename to use
% responses: cell array of TxMxR arrays, where T is number of time points, M is
%            dimensionality of features/states/etc., and R is the number of
%            simulated responses we have per feature set (i.e., if we
%            simulated responses for 3 different conditions, VNS, trauma,
%            and neutral, then R = 3). Each cell holds one participant's
%            responses
% Tx1 time vector to use to encode time

% (optional)
% 'vtn' or 'vt' or 'vn' or 'v' - If you input one of these strings, you
% are specifying that you want a legend on each plot with responses
% named according to input condition (e.g., 'vt' - VNS and trauma)

% 'verbose' if you want to see all warnings (default is not verbose)

% 'subjOverlay' - if this flag is provided, the next varargin cell needs to
% be an vector of subject numbers to use for a legend. What this flag
% entails is it will set meanResp = false. This will overlay all subjects'
% responses for a particular condition on top of each other while
% separating conditions into separate plots. If meanResp = true, means
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
    elseif strcmp(varargin{arg}, 'verbose')
        verbose = true;
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanResp = false;    % We will plot responses themselves
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
if ~exist('verbose', 'var')
    verbose = false;
end
if ~exist('meanResp', 'var')
    meanResp = true; % read documentation above for what this bool dictates
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
numSubjects = numel(responses);
numSamples = size(responses{1}, 1);
numFeats = size(responses{1}, 2);
numConds = size(responses{1}, 3);

% Group responses by condition; within each condition, stack by subject
resps_condGroup = cell(numConds, 1);
for cond_idx = 1:numConds
    % Initialize
    resps_condGroup{cond_idx} = zeros(numSamples, numFeats, numSubjects);
end

% Now fill the cell array's arrays
for cond_idx = 1:numConds
    for sub = 1:numSubjects
        resps_condGroup{cond_idx}(:, :, sub) = ...
            responses{sub}(:, :, cond_idx);
    end
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
    
    % Compute means within each condition (3rd dimension is subject)
    meanResps_condGroup = zeros(numSamples, numFeats, numConds);
    for cond_idx = 1:numConds
        meanResps_condGroup(:, :, cond_idx) = ...
            mean(resps_condGroup{cond_idx}, 3);
    end
    
    % Initialize stacked responses
    stackedResps = [];
    condLabels = [];    % To easily re-separate conditions
    
    % Loop and stack
    for cond_idx = 1:numConds
        stackedResps = [stackedResps; meanResps_condGroup(:, :, cond_idx)];
        condLabels = [condLabels; ...
            cond_idx*ones(size(meanResps_condGroup(:, :, cond_idx), 1), 1)];
    end
    
    % Check to make sure sizes are correct
    if size(stackedResps, 1) ~= numSamples*numConds || ...
            length(condLabels) ~= numSamples*numConds
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
    
    % Separate the PCA scores by subject and condition to make labeling of
    % plots easy (i.e., "unstack")
    plotReady = cell(numConds, 1);
    for cond_idx = 1:numConds
        plotReady{cond_idx} = ...
            pcaResult(condLabels == cond_idx, :);
        
        % Check size to make sure everything was stored properly
        if size(plotReady{cond_idx}, 1) ~= numSamples || ...
                size(plotReady{cond_idx}, 2) ~= numFeats
            error('plotReady sizing is incorrect')
        end
    end
    
    % Branch based on dimensions to plot
    if dim2plot == 2
        % Create figure
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Iterate through conditions and overlay
        hold on
        for cond_idx = 1:numConds
            % Plot one scatter point at a time to leverage time labels
            for t = 1:length(timeVec)
                s = scatter(plotReady{cond_idx}(t, 1), ...
                    plotReady{cond_idx}(t, 2), sizeLabels(t), ...
                    colorLabels(cond_idx, :), 'filled');
                s.MarkerFaceAlpha = alphaLabels(t);
                
                % If this is the final datapoint, then store legend
                % info
                if t == length(timeVec)
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
        legend(condNames, 'Location', 'best')
        legend('boxoff')
        hold off
        
        % Add title to entire figure
        title(plotTitle);
        
        % Save plot
        saveas(fig, fileName, 'fig');
        close all
        
    else
        % Create figure
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Iterate through conditions and overlay
        hold on
        for cond_idx = 1:numConds
            % Plot one scatter point at a time to leverage time labels
            for t = 1:length(timeVec)
                s = scatter3(plotReady{cond_idx}(t, 1), ...
                    plotReady{cond_idx}(t, 2), ...
                    plotReady{cond_idx}(t, 3), ...
                    sizeLabels(t), colorLabels(cond_idx, :), 'filled');
                s.MarkerFaceAlpha = alphaLabels(t);
                
                % If this is the final datapoint, then store legend
                % info
                if t == length(timeVec)
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
        
        % Add title to entire figure
        title(plotTitle);
        
        % Save plot
        saveas(fig, fileName, 'fig');
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
    % Create color map based on number of subjects
    rng(0);     % Set random seed so that colors are always the same
    colorLabels = rand(numSubjects, 3);
    
    % Create cell array of subject IDs based on subject numbers for legend
    subjectIDs = cell(length(subjectNums), 1);
    for sub = 1:numel(subjectIDs)
        subjectIDs{sub} = num2str(subjectNums(sub));
    end
    
    % Since I am not doing any averaging, I essentially need to stack and
    % flatten resps_condGroup to be a 2-D array of dimension
    % numSamples*numSubjects*numConds x numFeats. At the same time, I
    % will need to keep track of subjLabels and condLabels so I can later
    % separate when plotting
    subjLabels = [];
    condLabels = [];
    
    % Initialize the stacked result
    stackedResps = [];
    
    % Loop and stack
    for cond_idx = 1:numConds
        for sub_idx = 1:numSubjects
            stackedResps = [stackedResps; ...
                resps_condGroup{cond_idx}(:, :, sub_idx)];
            subjLabels = [subjLabels; subjectNums(sub_idx)*...
                ones(size(resps_condGroup{cond_idx}(:, :, sub_idx), 1), ...
                1)];
            condLabels = [condLabels; cond_idx*...
                ones(size(resps_condGroup{cond_idx}(:, :, sub_idx), 1), ...
                1)];
        end
    end
    
    % Check to make sure sizes are correct
    if size(stackedResps, 1) ~= numSamples*numSubjects*numConds || ...
            length(subjLabels) ~=  numSamples*numSubjects*numConds || ...
            length(condLabels) ~= numSamples*numSubjects*numConds
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
    
    % Separate the PCA scores by subject and condition to make labeling of
    % plots easy (i.e., "unstack")
    plotReady = cell(numConds, numSubjects);
    for cond_idx = 1:numConds
        for sub_idx = 1:numSubjects
            plotReady{cond_idx, sub_idx} = ...
                pcaResult(subjLabels == subjectNums(sub_idx) & ...
                condLabels == cond_idx, :);
            
            % Check size to make sure everything was stored properly
            if size(plotReady{cond_idx, sub_idx}, 1) ~= numSamples || ...
                    size(plotReady{cond_idx, sub_idx}, 2) ~= numFeats
                error('plotReady sizing is incorrect')
            end
        end
    end
    
    % Branch based on dimensions to plot
    if dim2plot == 2
        % Create figure
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Initialize tight subplot with a subplot per condition arranged
        % horizontally
        subs = tight_subplot(1, numConds);
        
        % Iterate through conditions
        for cond_idx = 1:numConds
            % Each condition will be on its own subplot
            axes(subs(cond_idx));
            hold on
            % Iterate through participants and overlay
            for sub_idx = 1:numSubjects
                % Plot one scatter point at a time to leverage time labels
                for t = 1:length(timeVec)
                    s = scatter(plotReady{cond_idx, sub_idx}(t, 1), ...
                        plotReady{cond_idx, sub_idx}(t, 2), sizeLabels(t), ...
                        colorLabels(sub_idx, :), 'filled');
                    s.MarkerFaceAlpha = alphaLabels(t);
                    
                    % If this is the final datapoint, then store legend
                    % info
                    if t == length(timeVec)
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
            if cond_idx == 3    % (I only want legend on rightmost plot)
                legend(subjectIDs, 'Location', 'eastoutside')
            end
            
            % Title it according to condition
            title(condNames{cond_idx});
            
            hold off
        end
        
        % Add title to entire figure
        sgtitle(plotTitle);
        
        % Save plot
        saveas(fig, fileName, 'fig');
        close all
        
    else
        % Create figure
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Initialize tight subplot with a subplot per condition arranged
        % horizontally
        subs = tight_subplot(1, numConds);
        
        % Iterate through conditions
        for cond_idx = 1:numConds
            % Each condition will be on its own subplot
            axes(subs(cond_idx));
            hold on
            % Iterate through participants and overlay
            for sub_idx = 1:numSubjects
                % Plot one scatter point at a time to leverage time labels
                for t = 1:length(timeVec)
                    s = scatter3(plotReady{cond_idx, sub_idx}(t, 1), ...
                        plotReady{cond_idx, sub_idx}(t, 2), ...
                        plotReady{cond_idx, sub_idx}(t, 3), ... 
                        sizeLabels(t), ...
                        colorLabels(sub_idx, :), 'filled');
                    s.MarkerFaceAlpha = alphaLabels(t);
                    
                    % If this is the final datapoint, then store legend
                    % info
                    if t == length(timeVec)
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
            
            if cond_idx == 3    % (I only want legend on rightmost plot)
                legend(subjectIDs, 'Location', 'eastoutside')
            end
            
            % Title it according to condition
            title(condNames{cond_idx});
            
            hold off
        end
        
        % Add title to entire figure
        sgtitle(plotTitle);
        
        % Save plot
        saveas(fig, fileName, 'fig');
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

