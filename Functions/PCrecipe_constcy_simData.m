function PCrecipe_constcy_simData(plotTitle, fileName, responses_cellSep, ...
    varargin)
% This function will take a set of simulated responses, apply PCA to each
% of them separately, and then aggregate the recipes produced to try and
% understand how consistent the first 3 PC recipes are

% ----------- Inputs ----------------- %
% plotTitle: first part of figure title to use for each figure
% fileName: first part of filename to use for each figure
% responses_cellSep: cell array of some data type, where principal
%                    component recipes will be computed for each cell and
%                    then aggregated to create consistency plots. Each cell
%                    can either contain a struct or another cell array. If
%                    another cell array, that just means we are dealing
%                    with an entire group's data, so the cell array within
%                    will be organized to contain a struct in each cell.
%                    Each struct will have R fields, where R represents the
%                    number of conditions (e.g., R = 3 for VNS, trauma, and
%                    neutral). Each field will contain a T_i x M array. T_i
%                    represents the number of timepoints for the ith
%                    condition (i \in {1, ..., R}), and M is the number of
%                    features

% (optional)
% 'verbose' - if you want to see all warnings (default is not verbose)

% 'subjOverlay' - if this flag is provided and responses_cellSep contains
%                 cell arrays within each cell, the PC recipes will be
%                 computed not on the means responses of each group, but
%                 rather on the concatenation of all subjects' responses

% 'fixPositive' - if this flag is provided, that means that we want to
%                  one of the feature directions either positive or
%                  negative for all the principal components to make them
%                  comparable in some way across cells. After this flag,
%                  the feature index that needs to remain positive should
%                  be provided

% 'featNames' - if this flag is provided, the next input should be a cell
%               array of feature names to be used as x-tick labels for
%               labeling x axis of consistency plots

% 'grpNames' - if this flag is provided, the next input should be a cell
%              array of subject/group names to be used as y-tick labels for
%              labeling y axis of consistency plots


% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'verbose')
        verbose = true; % We want the warnings
    elseif strcmp(varargin{arg}, 'subjOverlay')
        meanResp = false;   % We will consider each group's subjects separately
    elseif strcmp(varargin{arg}, 'fixPositive')
        featPositive = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'featNames')
        labelFeats = true;  % Set boolean to true
        featNames = varargin{arg + 1};  % Store feature names
    elseif strcmp(varargin{arg}, 'grpNames')
        labelGrps = true;   % Set boolean to true
        grpNames = varargin{arg + 1};   % Store group names
    end
end

% Set defaults
if ~exist('verbose', 'var'); verbose = false; end
if ~exist('meanResp', 'var'); meanResp = true; end
if ~exist('featPositive', 'var'); featPositive = -1; end
if ~exist('labelFeats', 'var'); labelFeats = false; end
if ~exist('labelGrps', 'var'); labelGrps = false; end

% We suppress these expected warnings if verbosity is not desired
if ~verbose
    % Specify the warning IDs to turn off
    wrnID = {'stats:pca:ColRankDefX'};
    
    % Turn these warnings off
    for id = 1:numel(wrnID)
        warning('off', wrnID{id})
    end
end

% Just to have number of cells (and subcells, if applicable) saved away
numCells = numel(responses_cellSep);
numSubcells = ones(numCells, 1);    % Set to one unless proven otherwise
for c = 1:numCells
    % If we're not dealing with subcells, we just keep # subcells = 1
    if iscell(responses_cellSep{c})
        numSubcells(c) = numel(responses_cellSep{c});
    end
end

% Check to see how many conditions and features; store away condition names
% while I'm at it
if numSubcells(1) > 1
    condNames = fieldnames(responses_cellSep{1}{1});
    numConds = length(condNames);
    numFeats = size(responses_cellSep{1}{1}.(condNames{1}), 2);
else
    condNames = fieldnames(responses_cellSep{1});
    numConds = length(condNames);
    numFeats = size(responses_cellSep{1}.(condNames{1}), 2);
end

% Store away number of samples in each condition's response so we know how
% to initilaize arrays later
numSamples = zeros(numConds, 1);
for cond = 1:numConds
    if numSubcells(1) > 1
        numSamples(cond) = ...
            size(responses_cellSep{1}{1}.(condNames{cond}), 1);
    else
        numSamples(cond) = ...
            size(responses_cellSep{1}.(condNames{cond}), 1);
    end
end



% For each cell, I now need to form a stacked matrix across all timepoints,
% conditions, and subjects if applicable
stackedArrays = cell(numCells, 1);      % Initialize cell array for result
for c = 1:numCells
    if numSubcells(c) > 1
        % This means we are working with multiple subjects' data
        
        % First, let's separate responss by condition; within
        % each condition, stack subjects in 3rd dimension
        cellSpec_resps_condSep = cell(numConds, 1);
        for cond = 1:numConds
            % Loop through subcells to stack subjects' data for this cond.
            for subcell = 1:numSubcells(c)
                cellSpec_resps_condSep{cond} = cat(3, ...
                    cellSpec_resps_condSep{cond}, ...
                    responses_cellSep{c}{subcell}.(condNames{cond}));
            end
        end
        
        % Initialize a new cell array to store intermediate result after
        % next step
        cellSpec_resps_condReady = cell(numConds, 1);
        % Branch based on subject overlaying or not
        if meanResp
            % In this case, we need to compute the mean condition-specific
            % responses and then concatenate
            for cond = 1:numConds
                cellSpec_resps_condReady{cond} = ...
                    mean(cellSpec_resps_condSep{cond}, 3);
            end
        else
            % In this case, we need to unstack along 3rd dimension and
            % instead concatenate along 1st dimension
            for cond = 1:numConds
                for indx = 1:size(cellSpec_resps_condSep{cond}, 3)
                    cellSpec_resps_condReady{cond} = ...
                        [cellSpec_resps_condReady{cond}; ...
                        cellSpec_resps_condSep{cond}(:, :, indx)];
                end
            end
        end
        
        % We can now stack across conditions and store away
        for cond = 1:numConds
            stackedArrays{c} = [stackedArrays{c}; ...
                cellSpec_resps_condReady{cond}];
        end
    else
        % In this case, I don't need to worry about stacking across
        % subcells. I only need to stack across conditions
        for cond = 1:numConds
            stackedArrays{c} = [stackedArrays{c}; ...
                responses_cellSep{c}.(condNames{cond})];
        end
    end
end

% Now, we can take the stackedArrays, apply PCA to each of them, and store
% away the first 3 PC recipes
PCrecipes = zeros(numCells, numFeats, 3);
varExplains = zeros(numCells, 3);
for c = 1:numCells
    pcaReady = zeros(size(stackedArrays{c}));
    
    % Unit variance
    for col = 1:size(stackedArrays{c}, 2)
        pcaReady(:, col) = ...
            (1/std(stackedArrays{c}(:, col)))*stackedArrays{c}(:, col);
    end
    
    % Apply PCA and store away PC recipes and variance explained
    [pcaLoads, ~, ~, ~, varExplain, ~] = pca(pcaReady, 'Centered', false);
    
    % Store recipes
    for dim = 1:3
        PCrecipes(c, :, dim) = pcaLoads(:, dim);
        varExplains(c, dim) = varExplain(dim);
    end
end

% If I need to fix a particular feature direction as positive across all PC
% recipes, we can do that now
if featPositive > -1
    for c = 1:numCells
        for dim = 1:3
            if PCrecipes(c, featPositive, dim) < 0
                PCrecipes(c, :, dim) = -1*PCrecipes(c, :, dim);
            end
        end
    end
end


% Also create a binarized version of the PC recipes that just stores sign
signRecipes = sign(PCrecipes);

% Take this binarized version and convert it into a form that tells you
% percentage up/down
percentRecipes = zeros(size(signRecipes, 2), size(signRecipes, 3));
for f = 1:numFeats
    for dim = 1:3
        vec = signRecipes(:, f, dim);
        if sum(vec) > 0
            % I need to count how many of the elements are +1 and then
            % divide by total elements to get a percent up
            countPos1 = length(vec(vec == 1));
            
            percentRecipes(f, dim) = 100*countPos1/length(vec);
        else
            % I need to count how many of the elements are -1 and then
            % divide by total elements to get a percent down
            countNeg1 = length(vec(vec == -1));
            
            percentRecipes(f, dim) = -100*countNeg1/length(vec);
        end
    end
end


% Compute mean and SEM for recipes because we will need later
meanRecipes = squeeze(mean(PCrecipes, 1));
semRecipes = squeeze((1/sqrt(size(PCrecipes, 1)))*std(PCrecipes, 0, 1));


% Now we can create the plots we would like to
for dim = 1:3
    % Create image figure
    fig = figure(1);
    set(fig, 'Visible', 'on');
    
    % Create image
    imagesc(PCrecipes(:, :, dim));
    
    % Use red blue colormap that goes from blue to white to red to make
    % directionality easy to see
    colormap(fig, redblue(256));
    
    % Add colorbar
    colorbar;
    
    % Add labels
    if labelFeats
        xlabel('Feature Name');
        xticks(1:numFeats);
        xticklabels(featNames);
        xtickangle(45);
    else
        xlabel('Feature Index');
    end
    
    % Add labels
    if labelGrps
        ylabel('ID')
        yticks(1:numCells);
        yticklabels(grpNames);
    else
        ylabel('Grouping Index')
    end
    
    
    % Add title
    title([plotTitle, ' - Recipes for PC ', num2str(dim)]);
    
    % Save plot
    saveas(fig, [fileName, '_RecipesImg_dim', num2str(dim)], 'fig');
    
    
    % Create mean +- 95% CI figure
    fig = figure(2);
    set(fig, 'Visible', 'on');
    
    % Create x coordinates for plotting purposes (will help with error bars
    % especially)
    xCoords = 1:numFeats;
    
    hold on;
    % Add bars for means
    bar(xCoords, meanRecipes(:, dim))
    
    % Add error bars for 95% CIs
    errorbar(xCoords, meanRecipes(:, dim), 1.96*semRecipes(:, dim), ...
        'LineStyle', 'none')
    hold off;
    
    % Add labels
    if labelFeats
        xlabel('Feature Name');
        xticks(1:numFeats);
        xticklabels(featNames);
        xtickangle(45);
    else
        xlabel('Feature Index');
    end
    
    % Add label
    ylabel('Coefficient')
    
    % Add title
    title([plotTitle, ' - Mean \pm 95% CI for PC ', num2str(dim)]);
    
    % Save plot
    saveas(fig, [fileName, '_Mean95CI_dim', num2str(dim)], 'fig');
    
    
    
    % Create an image using signs
    fig = figure(3);
    set(fig, 'Visible', 'on');
    
    % Create image
    imagesc(signRecipes(:, :, dim));
    
    % Use red blue colormap that goes from blue to white to red to make
    % directionality easy to see
    colormap(fig, redblue(3));
    
    % Add colorbar
    colorbar;
    
    % Add labels
    if labelFeats
        xlabel('Feature Name');
        xticks(1:numFeats);
        xticklabels(featNames);
        xtickangle(45);
    else
        xlabel('Feature Index');
    end
    
    % Add labels
    if labelGrps
        ylabel('ID')
        yticks(1:numCells);
        yticklabels(grpNames);
    else
        ylabel('Grouping Index')
    end
    
    
    % Add title
    title([plotTitle, ' - Recipe Signs ', num2str(dim)]);
    
    % Save plot
    saveas(fig, [fileName, '_RecipeSignsImg_dim', num2str(dim)], 'fig');
    
    
    
    % Create bar plot using percent up/down
    fig = figure(4);
    set(fig, 'Visible', 'on');
    
    % Create x coordinates for plotting purposes
    xCoords = 1:numFeats;
    
    % Add bars for percentages
    hold on;
    bar(xCoords, percentRecipes(:, dim))
    yline(-50, '--k');
    yline(50, '--k');
    hold off;
    
    % Add labels
    if labelFeats
        xlabel('Feature Name');
        xticks(1:numFeats);
        xticklabels(featNames);
        xtickangle(45);
    else
        xlabel('Feature Index');
    end
    
    % Add label
    ylabel('Percentage Positive or Negative')
    
    % Add title
    title([plotTitle, ' - Percentage Up/Down for PC ', num2str(dim)]);
    
    % Save plot
    saveas(fig, [fileName, '_PercentPosNeg_dim', num2str(dim)], 'fig');
    
    
    
    % Create histogram of variance explained
    fig = figure(5);
    set(fig, 'Visible', 'on');
    
    % Create histogram for this dimension
    histogram(varExplains(:, dim))
    
    % Label axes
    xlabel('Variance Explained (%)')
    ylabel('Frequency')
    
    % Add title
    title([plotTitle, ' - Variance Explained for PC ', num2str(dim)]);
    
    % Save plot
    saveas(fig, [fileName, '_VarExplain_dim', num2str(dim)], 'fig');
    
    
    % Close figures
    close all
end

end

