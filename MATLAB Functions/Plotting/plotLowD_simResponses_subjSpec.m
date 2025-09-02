function plotLowD_simResponses_subjSpec(dim2plot, plotTitle, fileName, ...
    responses, timeVec, varargin)
% This function will take simulated responses, apply PCA to them, and then
% plot and save them according to inputs provided

% ---------- Inputs -------------------- %
% dim2plot: scalar = 2 or 3 that dictates whether to plot in 2-D or 3-D
% plotTitle: figure title to use
% fileName: filename to use
% responses: TxMxR array, where T is number of time points, M is
%            dimensionality of features/states/etc., and R is the number of
%            simulated responses we have per feature set (i.e., if we
%            simulated responses for 3 different conditions, VNS, trauma,
%            and neutral, then R = 3)
% Tx1 time vector to use to encode time

% (optional)
% 'vtn' or 'vt' or 'vn' or 'v' - If you input one of these strings, you
% are specifying that you want a legend on each plot with responses
% named according to input condition (e.g., 'vt' - VNS and trauma)

% 'verbose' if you want to see all warnings (default is not verbose)

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
    end
end

% Set defaults
if ~exist('condNames', 'var')
    % If they don't provide condition names, I'll just number them
    condNames = {};
    for r = 1:size(responses, 3)
        condNames = [condNames; {['Cond. ', num2str(r)]}]; 
    end
end
if ~exist('verbose', 'var')
    verbose = false;
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
numConds = size(responses, 3);
numFeats = size(responses, 2);
numSamples = size(responses, 1);

% Stack condition responses on top of each other, if needed
if numConds > 1
    % Initialize
    stackedResps = zeros(numSamples*numConds, numFeats);
    
    % Loop through and stack
    for cond_idx = 1:numConds
        stackedResps(((cond_idx-1)*numSamples + 1):cond_idx*numSamples, :) = ...
            responses(:, :, cond_idx);
    end
    
else
    % Store as is
    stackedResps = responses;
end

% Each column should be unit variance
% (I do not make zero mean because the simulated responses start at the
% origin since they're initialized there and I want to keep starting point
% at origin)
pcaReady = zeros(size(stackedResps));    % Initialize
for col = 1:size(stackedResps, 2)
    pcaReady(:, col) = (1/std(stackedResps(:, col)))*stackedResps(:, col);
end

% Apply PCA and store away variance explained and PCs
% (outputs I ignore are the principal component coefficients, the variances
% in the latent space, Hotelling's T-squared statistic for each
% observation, and the mean of each column)
[pcaLoads, pcaResult, ~, ~, varExplain, ~] = pca(pcaReady, ...
            'Centered',false);

% Separate the responses back to the way they were input (i.e., "unstack")
plotReady = zeros(size(responses));
for cond_idx = 1:numConds
    plotReady(:, :, cond_idx) = ...
        pcaResult(((cond_idx-1)*numSamples + 1):cond_idx*numSamples, :);
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

% Create basic color labels (blue, red, green)
colorLabels = (1/255)*[221, 85, 255; ...
    255, 85, 85; ...
    95, 211, 188];

if dim2plot == 2
    % 2-D scatter plot
    fig = figure(1);
    set(fig, 'Visible', 'on');
    
    hold on
    % Iterate through conditions
    for cond_idx = 1:numConds
        for t = 1:length(timeVec)
            % Plot one scatter point at a time
            s = scatter(plotReady(t, 1, cond_idx), ... 
                plotReady(t, 2, cond_idx), sizeLabels(t), ...
                colorLabels(cond_idx, :), 'filled');
            s.MarkerFaceAlpha = alphaLabels(t);
            
            % If this is the final datapoint, then store legend info
            if t == length(timeVec)
                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
            else
                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
        end
    end
    hold off
    
    % Include variance explained by each of the principal components in the
    % axis labels
    xlabel(['PC1 Variance Explained : ', num2str(varExplain(1)), '%']);
    ylabel(['PC2 Variance Explained : ', num2str(varExplain(2)), '%']);
    
    % Title it according to input provided
    title(plotTitle);
    
    % Add a legend
    legend(condNames, 'Location', 'best');
    legend('boxoff')
    
    savefig(fig, fileName);
    close all
else
    % 3-D scatter plot
    fig = figure(1);
    set(fig, 'Visible', 'on');
    
    hold on
    % Iterate through conditions
    for cond_idx = 1:numConds
        for t = 1:length(timeVec)
            % Plot one scatter point at a time
            s = scatter3(plotReady(t, 1, cond_idx), ...
                plotReady(t, 2, cond_idx), plotReady(t, 3, cond_idx), ...
                sizeLabels(t), colorLabels(cond_idx, :), 'filled');
            s.MarkerFaceAlpha = alphaLabels(t);
            
            % If this is the final datapoint, then store legend info
            if t == length(timeVec)
                s.Annotation.LegendInformation.IconDisplayStyle = 'on';
            else
                s.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            
        end
    end
    hold off
    
    % Include variance explained by each of the principal components in the
    % axis labels
    xlabel(['PC1 Variance Explained : ', num2str(varExplain(1)), '%']);
    ylabel(['PC2 Variance Explained : ', num2str(varExplain(2)), '%']);
    zlabel(['PC3 Variance Explained : ', num2str(varExplain(3)), '%']);
    
    % Title it according to input provided
    title(plotTitle);
    
    % Add a legend
    legend(condNames, 'Location', 'best');
    legend('boxoff')
    
    savefig(fig, fileName);
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

