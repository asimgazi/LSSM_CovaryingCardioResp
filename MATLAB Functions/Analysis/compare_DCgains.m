function compare_DCgains(plotTitle, fileName, Aparams1, Aparams2, ...
    Cparams1, Cparams2, featNames, Bparams1, Bparams2, inputsOpt, groupNames)
% This function is passed parameters from a C matrix (state-to-output)
% terms and terms from a B matrix (input-to-state terms)
% for two groups' models and compares the C*B products, feature by feature,
% between the subjects in both groups (unpaired t-tests or ranksum tests)

% Note: I didn't bother changing too many of the comments and names from the
% compare_CBproducts function, so many of the variable names and comments
% may not match up

% ---------- Inputs ---------- %
% plotTitle: Just so the plots have some meaningful title
% fileName: file name to use when saving figure
% Aparams1: a cell array of A matrices for group 1, one cell per subject
% Aparams2: a cell array of A matrices for group 2, one cell per subject
% Cparams1: a cell array of C matrices for group 1, one cell per subject
% Cparams2: a cell array of C matrices for group 2, one cell per subject
% featNames: Cell array of feature names to use when labeling
% Bparams1: a cell array of B matrices for group 1, one cell per subject
% Bparams2: a cell array of B matrices for group 2, one cell per subject
% inputsOpt: string that specifies which protocol conditions
%            inputs were created for -
%            'v': just VNS
%            'vt': VNS and trauma recall
%            'vn': VNS and neutral condition
%            'vtn': VNS, trauma recall, and neutral
% groupNames: 1x2 cell array storing each of the two groups' names 

% Parse varargin
if strcmp(inputsOpt, 'vtn')
    numInputs = 3;  % Store number of inputs
elseif strcmp(inputsOpt, 'vt') || strcmp(inputsOpt, 'vn')
    numInputs = 2;  % Store number of inputs
else
    numInputs = 1;  % Store number of inputs
end

% We arbitrarily use the first group's first subject's C matrix to check
% how many features
numFeats = size(Cparams1{1}, 1);

% Initialize arrays to store the products for both groups
% I am going to transpose the products to have each input as its own
% row, each feature as its own column
% I am going to stack the subjects' inputs in 3rd dimension
CBproducts1 = zeros(numInputs, numFeats, numel(Cparams1));
CBproducts2 = zeros(numInputs, numFeats, numel(Cparams2));

% Compute products and stack
for sub = 1:numel(Cparams1)
    CBproducts1(:, :, sub) = ...
        transpose(Cparams1{sub}*...
        inv(eye(size(Aparams1{sub}, 1)) - Aparams1{sub})*Bparams1{sub});
end
for sub = 1:numel(Cparams2)
    CBproducts2(:, :, sub) = ...
        transpose(Cparams2{sub}*...
        inv(eye(size(Aparams2{sub}, 1)) - Aparams2{sub})*Bparams2{sub});
end

% Now, I can create a set of boxplots for each feature, one boxplot per
% input
for f = 1:numFeats
    % Create an overall figure for this feature
    fig = figure(1);
    set(fig, 'Visible', 'on');
    for i = 1:numInputs
        % Access this subplot for this input
        subplot(1, numInputs, i);
        
        % Let's set aside our two vectors for easy access
        vec1 = squeeze(CBproducts1(i, f, :)); 
        vec2 = squeeze(CBproducts2(i, f, :));
        
        % Let's compute a p-value using the appropriate test
        if swtest(vec1) || swtest(vec2)
            pVal = ranksum(vec1, vec2);
        else
            [~, pVal] = ttest2(vec1, vec2);
        end
        
        % To make boxplotting and plotting raw data easier, I am going to
        % reoorgnize things
        boxPlotReady = zeros(length(vec1) + length(vec2), 2);
        boxPlotReady(:, 1) = [vec1; vec2];
        boxPlotReady(:, 2) = [ones(length(vec1), 1); 2*ones(length(vec2), 1)];
        
        % Create boxplot
        hold on;
        boxplot(boxPlotReady(:, 1), boxPlotReady(:, 2));
        scatter(boxPlotReady(:, 2), boxPlotReady(:, 1), 'o');
        hold off
        
        % Label y axis
        ylabel('DC Gain Value')
        
        % Label x axis
        xlabel('Group Names');
        xticks([1, 2]);
        xticklabels(groupNames);
        
        % Title the subplot
        if i == 1
            title(['VNS/Sham Term; P = ', num2str(pVal)])
        elseif i == 2
            if strcmp(inputsOpt, 'vtn') || strcmp(inputsOpt, 'vt')
                title(['Trauma Term; P = ', num2str(pVal)])
            else
                title(['Neutral Term; P = ', num2str(pVal)])
            end
        else
            title(['Neutral Term; P = ', num2str(pVal)])
        end
    end
    % Title the entire plot and save
    sgtitle([plotTitle, ' - ', featNames{f}])
    savefig(fig, [fileName, '_', featNames{f}]);
    
    close all
end


end

