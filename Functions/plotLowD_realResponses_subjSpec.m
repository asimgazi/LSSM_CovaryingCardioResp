function plotLowD_realResponses_subjSpec(dim2plot, subjectID, saveDir, ...
    featSegm, Fs, varargin)
% This function will take real responses, apply PCA to them, and then
% plot and save them according to inputs provided

% ---------- Inputs -------------------- %
% dim2plot: scalar = 2 or 3 that dictates whether to plot in 2-D or 3-D
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
%               'condSep'
%               'condAvg'
%               'condAvg_train': if this option is selected, then the next
%               inputs should specify the train-test
%               split condition code (e.g., 6.5 for split between 6 and 7)

% 'verbose' if you want to see all warnings (default is not verbose)

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
condNames = fieldnames(featSegm);
numConds = length(condNames);

if condAvg
    % Initialize a condition-averaged response struct to store result
    feat_condAvg = struct();
    for cond = 1:numConds
        % Initialized to NaNs as default in case no data are available to
        % replace (will ignore NaN responses during plotting)
        feat_condAvg.(condNames{cond}) = ...
            nan*ones(size(featSegm.(condNames{cond}), 1), ...
            size(featSegm.(condNames{cond}), 2));
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
    
    % Now that we have the condition averages for this subject, we can
    % stack the condition responses on top of each other, if needed
    stackedResps = [];
    legendNames = {};   % Keeps track of conditions that are not NaN
    condSamples = [];   % Keeps track of # samples per condition
    % Loop through and stack if not NaN
    for cond = 1:numConds
        if ~isnan(sum(feat_condAvg.(condNames{cond}), 'all'))
            stackedResps = [stackedResps; feat_condAvg.(condNames{cond})];
            
            % Store this condition name away
            legendNames{cond} = condNames{cond};
            
            % Store the number of samples away
            condSamples = [condSamples, ...
                size(feat_condAvg.(condNames{cond}), 1)];
        end
    end
    
    % Each column should be unit variance
    % (I do not make zero mean because I am interested in biases from
    % equilibrium)
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
    % I use condSamples and legendNames because they skipped missing conditions
    plotReady = struct();
    for cond = 1:length(condSamples)
        plotReady.(legendNames{cond}) = ...
            pcaResult((sum(condSamples(1:(cond-1))) + 1):...
            (sum(condSamples(1:(cond-1))) + condSamples(cond)), :);
    end
    
    % Create basic color labels (blue, red, green) for protocol conditions
    colorLabels.v = (1/255)*[221, 85, 255];
    colorLabels.t = (1/255)*[255, 85, 85];
    colorLabels.n = (1/255)*[95, 211, 188];
    
    % To encode time, I am going to create a time vector according to the
    % longest time protocol condition we have and then just use a subset of
    % that for the other conditions
    timeVec = 0:(1/Fs):(max(condSamples)-1)/Fs;
    
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
    
    if dim2plot == 2
        % 2-D scatter plot
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        hold on
        % Iterate through conditions
        for cond = 1:length(condSamples)
            for t = 1:size(plotReady.(legendNames{cond}), 1)
                % Plot one scatter point at a time
                s = scatter(plotReady.(legendNames{cond})(t, 1), ...
                    plotReady.(legendNames{cond})(t, 2), sizeLabels(t), ...
                    colorLabels.(legendNames{cond}), 'filled');
                s.MarkerFaceAlpha = alphaLabels(t);
                
                % If this is the final datapoint, then store legend info
                if t == size(plotReady.(legendNames{cond}), 1)
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
        title(subjectID);
        
        % Add a legend
        legend(legendNames, 'Location', 'best');
        legend('boxoff')
        
        savefig(fig, [saveDir, 'Subject ', subjectID]);
        close all
    else
        % 3-D scatter plot
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        hold on
        % Iterate through conditions
        for cond = 1:length(condSamples)
            for t = 1:size(plotReady.(legendNames{cond}), 1)
                % Plot one scatter point at a time
                s = scatter3(plotReady.(legendNames{cond})(t, 1), ...
                    plotReady.(legendNames{cond})(t, 2), ...
                    plotReady.(legendNames{cond})(t, 3), sizeLabels(t), ...
                    colorLabels.(legendNames{cond}), 'filled');
                s.MarkerFaceAlpha = alphaLabels(t);
                
                % If this is the final datapoint, then store legend info
                if t == size(plotReady.(legendNames{cond}), 1)
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
        title(subjectID);
        
        % Add a legend
        legend(legendNames, 'Location', 'best');
        legend('boxoff')
        
        savefig(fig, [saveDir, 'Subject ', subjectID]);
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
    saveas(fig, [saveDir, 'Subject ', subjectID, '_PCloadings'], 'fig');
    
    % Close all plots
    close all
    
else
    % In this case, we do not want to average across repetitions of a
    % condition, so each available repetition will be treated on its own
    
    % This function is created for the DARPA PTSD dataset, so the
    % legend naming will be in accordance with the protocol ordering:
    % |N1|N2|T1V1|T2V2|V3|V4|N3|N4|T3V5|T4V6|
    % (10 total conditions; 14 total responses)
    % Defining the legend names in this way will help us color code all
    % participants' responses equivalently, taking into consideration
    % missing responses
    legendNames.v = {'Stim 1', 'Stim 2', 'Stim 3', 'Stim 4', 'Stim 5', 'Stim 6'};
    legendNames.t = {'Trauma 1', 'Trauma 2', 'Trauma 3', 'Trauma 4'};
    legendNames.n = {'Neutral 1', 'Neutral 2', 'Neutral 3', 'Neutral 4'};
    
    % Specify color labels in a similar way to legend names
    colorLabels.v = (1/255)*...
        [141, 95, 211; ...
        212, 42, 255; ...
        255, 85, 221; ...
        68, 33, 120; ...
        222, 135, 205; ...
        85, 0, 212];
    colorLabels.t = (1/255)*...
        [255, 42, 127; ...
        255, 85, 85; ...
        211, 95, 95; ...
        211, 95, 141];
    colorLabels.n = (1/255)*...
        [85, 212, 0; ...
        135, 222, 135; ...
        44, 160, 90; ...
        95, 211, 188];
    
    % Since we do not need to average anything, we can directly start with
    % stacking the condition responses on top of each other, if needed
    stackedResps = [];
    condSamples = zeros(1, numConds);   % Keeps track of # samples per condition
    % Loop through and stack if not NaN
    for cond = 1:numConds
        for rep = 1:size(featSegm.(condNames{cond}), 3)
            if ~isnan(sum(featSegm.(condNames{cond})(:, :, rep), 'all'))
                % Stack
                stackedResps = [stackedResps; featSegm.(condNames{cond})(:, :, rep)];
            end
        end
        % Loop backwards so indexing works
        for rep = size(featSegm.(condNames{cond}), 3):-1:1
            if isnan(sum(featSegm.(condNames{cond})(:, :, rep), 'all'))
                % Remove this repetition from legend and color labels
                legendNames.(condNames{cond})(rep) = [];
                colorLabels.(condNames{cond})(rep, :) = [];
            end
        end
        % Store the number of samples for this condition so we know how to
        % unstack once dim reduction is done
        condSamples(cond) = size(featSegm.(condNames{cond}), 1);
    end
    
    % Each column should be unit variance
    % (I do not make zero mean because I am interested in biases from
    % equilibrium)
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
    % I use legendNames and condSamples to my advantage
    plotReady = struct();
    counter = 1;
    for cond = 1:numConds
        for rep = 1:length(legendNames.(condNames{cond}))
            plotReady.(condNames{cond})(:, :, rep) = ...
                pcaResult(counter:(counter - 1 + condSamples(cond)), :);
            
            % Increment counter
            counter = counter + condSamples(cond);
        end
    end
    
    % To encode time, I am going to create a time vector according to the
    % longest time protocol condition we have and then just use a subset of
    % that for the other conditions
    timeVec = 0:(1/Fs):(max(condSamples)-1)/Fs;
    
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
    
    if dim2plot == 2
        % 2-D scatter plot
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        hold on
        % Iterate through conditions
        for cond = 1:numConds
            % Iterate through available repetitions
            for rep = 1:length(legendNames.(condNames{cond}))
                
                for t = 1:size(plotReady.(condNames{cond}), 1)
                    % Plot one scatter point at a time
                    s = scatter(plotReady.(condNames{cond})(t, 1, rep), ...
                        plotReady.(condNames{cond})(t, 2, rep), sizeLabels(t), ...
                        colorLabels.(condNames{cond})(rep, :), 'filled');
                    s.MarkerFaceAlpha = alphaLabels(t);
                    
                    % If this is the final datapoint, then store legend info
                    if t == size(plotReady.(condNames{cond}), 1)
                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    else
                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                end
                
            end
        end
        hold off
        
        % Include variance explained by each of the principal components in the
        % axis labels
        xlabel(['PC1 Variance Explained : ', num2str(varExplain(1)), '%']);
        ylabel(['PC2 Variance Explained : ', num2str(varExplain(2)), '%']);
        
        % Title it according to input provided
        title(subjectID);
        
        % Add a legend
        concatLegend = {};
        for cond = 1:numConds
            % Create new cell array with available protocol info by
            % concatenating the remaining legend names
            concatLegend = [concatLegend, legendNames.(condNames{cond})];
        end
        legend(concatLegend, 'Location', 'best');
        legend('boxoff')
        
        savefig(fig, [saveDir, 'Subject ', subjectID]);
        close all
    else
        % 3-D scatter plot
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        hold on
        % Iterate through conditions
        for cond = 1:numConds
            % Iterate through available repetitions
            for rep = 1:length(legendNames.(condNames{cond}))
                
                for t = 1:size(plotReady.(condNames{cond}), 1)
                    % Plot one scatter point at a time
                    s = scatter3(plotReady.(condNames{cond})(t, 1, rep), ...
                        plotReady.(condNames{cond})(t, 2, rep), ...
                        plotReady.(condNames{cond})(t, 3, rep), ...
                        sizeLabels(t), ...
                        colorLabels.(condNames{cond})(rep, :), 'filled');
                    s.MarkerFaceAlpha = alphaLabels(t);
                    
                    % If this is the final datapoint, then store legend info
                    if t == size(plotReady.(condNames{cond}), 1)
                        s.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    else
                        s.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
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
        title(subjectID);
        
        % Add a legend
        concatLegend = {};
        for cond = 1:numConds
            % Create new cell array with available protocol info by
            % concatenating the remaining legend names
            concatLegend = [concatLegend, legendNames.(condNames{cond})];
        end
        legend(concatLegend, 'Location', 'best');
        legend('boxoff')
        
        savefig(fig, [saveDir, 'Subject ', subjectID]);
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
    saveas(fig, [saveDir, 'Subject ', subjectID, '_PCloadings'], 'fig');
    
    % Close all plots
    close all
    
end

end

