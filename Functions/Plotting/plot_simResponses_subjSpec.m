function plot_simResponses_subjSpec(featPer, titlePrefix, fileNamePrefix, ...
    responses, timeVec, varargin)
% This function takes simulated responses and plots and saves them
% according to inputs provided

% --------- Inputs --------------- %
% featPer: number of feature subplots per figure
% titlePrefix: figure title prefix to use; the suffix will be _1, _2, and 
%              so on based on (# features/featPer)
% fileNamePrefix: filename prefix to use; the suffix will be _1, _2, and so
%                 on based on (# features/featPer)
% responses: TxMxR array, where T is number of time points, M is
%            dimensionality of features/states/etc., and R is the number of
%            simulated responses we have per feature set (i.e., if we
%            simulated responses for 3 different conditions, VNS, trauma,
%            and neutral, then R = 3)
% timeVec: Tx1 time vector to plot against (in seconds)

% (optional)
% 'vtn' or 'vt' or 'vn' or 'v' - If you input one of these strings, you
% are specifying that you want a legend on each plot based on the overlaid
% responses and that the responses will be named according to input cond.
% (e.g., 'vt' - VNS and trauma)

% 'names' - if this flag is provided, the next varargin cell needs to be a
% cell array of names corresponding to the simulated responses provided. If
% no names are provided, the responses are simply labeled by index

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
if ~exist('responseNames', 'var')
    % If they don't provide response names, I'll just number them
    responseNames = {};
    for m = 1:size(responses, 2)
        responseNames = [responseNames; {['Dim. ', num2str(m)]}]; 
    end
end

% Just so we have them as easy to understand variables
numFeats = size(responses, 2);
numConds = size(responses, 3);

% Find out how many features we're dealing with and figures we'll need
numFigs = ceil(numFeats/featPer);

% Start creating figures
feat_idx = 1;           % To help us iterate through features/states/etc.
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
            % Iterate through conditions and overlay
            hold on
            for cond_idx = 1:numConds
                plot(timeVec, responses(:, feat_idx, cond_idx))
                ylabel(responseNames{feat_idx});
                xlabel('Time (s)');
            end
            
            % Add legend if needed
            if numel(condNames) > 1
                legend(condNames);
            end
            hold off
        end
        
        % Increment feature index
        feat_idx = feat_idx + 1;
    end
    
    % Add title to entire subplot grid
    sgtitle([titlePrefix, '_', num2str(fig_idx)]);
    
    % Save plot
    saveas(fig, [fileNamePrefix, '_', num2str(fig_idx)], ...
        'fig');
    close all
end


end

