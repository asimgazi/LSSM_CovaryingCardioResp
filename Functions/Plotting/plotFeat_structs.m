function plotFeat_structs(featPer, saveDir, varargin)
% This function takes structs and plots them according to the options
% provided

% ------ Input Instructions ------- %
% featPer: features per plot
% saveDir: the directory to save the plots to
% 1. One of the inputs needs to be a string that is either 'featOverlay' or
%    'featSeparate'. This string will dictate whether structs of features
%    are overlaid on top of each other or kept in separate figures. By
%    default, each feature within the struct will be on a separate subplot,
%    but this option is more for a new set of features and an old set of
%    features, for example. This can be used to understand the effects of a
%    particular operation on the features
% 2. If you want to input a struct of features, you need to either specify
%    whether the set of features is 'featNew' or 'featOld'. This function,
%    as of 1/8/2022, allows for two sets of features max, where these
%    features can either be plotted separately or overlaid. To do this, you
%    first input the string 'featNew' and then you input the features next.

% -------- Use examples -------- %
% plotFeat_structs(6, ['Test Figures', filesep], 'featOverlay', 'featNew', newFeats, ...
% 'featOld', oldFeats);
% plotFeat_structs(12, ['Wow Figures', filesep], 'featSeparate', 'featNew', newFeats);


% We first need to parse the varargin looking for whether to overlay
% features if provided with more than one struct
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'featOverlay')
        overlayFeats = true;
    elseif strcmp(varargin{i}, 'featNew')
        newFeats = varargin{i+1};
    elseif strcmp(varargin{i}, 'featOld')
        oldFeats = varargin{i+1};
    end
end

% Set defaults and booleans to know how to proceed later
if ~exist('overlayFeats', 'var')
    overlayFeats = false;
end
if overlayFeats && exist('newFeats', 'var') && exist('oldFeats', 'var')
    twoPlots = false;
elseif exist('newFeats', 'var') && exist('oldFeats', 'var')
    twoPlots = true;
else
    twoPlots = false;
    
    % Assign the single structure of features so I can easily work with it
    % later without checking for which one exists
    if exist('newFeats', 'var')
        singleStruct = newFeats;
        singleNames = fieldnames(newFeats);
    else
        singleStruct = oldFeats;
        singleNames = fieldnames(oldFeats);
    end
end


% We now know how to proceed! :)


% Split procedures based on two plots or one
if twoPlots
    % Save feature names for plot labeling
    newNames = fieldnames(newFeats);
    oldNames = fieldnames(oldFeats);
    
    % Find how many features we're dealing with
    newNum = length(newNames); newPlots = ceil(newNum/featPer);
    oldNum = length(oldNames); oldPlots = ceil(oldNum/featPer);
    
    % Take care of new
    feat_idx = 1;
    for plot_idx = 1:newPlots
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Initialize a tight subplot
        subs = tight_subplot(featPer, 1);
        
        for sub_idx = 1:featPer
            if feat_idx <= newNum
                % Plot in subplot
                axes(subs(sub_idx));
                plot(newFeats.(newNames{feat_idx})(:, 1), ...
                    newFeats.(newNames{feat_idx})(:, 2), 'k.')
                ylabel(newNames{feat_idx});
                xlabel('Time (s)');
            end
            
            % Add to feat_idx
            feat_idx = feat_idx + 1;
        end
        
        % Save plot
        saveas(fig, [saveDir, 'New Feat Figure ', num2str(plot_idx)], ...
            'fig');
        close all
    end
    
    % Take care of old
    feat_idx = 1;
    for plot_idx = 1:oldPlots
        fig = figure(1);
        set(fig, 'Visible', 'on');
        
        % Initialize a tight subplot
        subs = tight_subplot(featPer, 1);
        
        for sub_idx = 1:featPer
            if feat_idx <= oldNum
                % Plot in subplot
                axes(subs(sub_idx));
                plot(oldFeats.(oldNames{feat_idx})(:, 1), ...
                    oldFeats.(oldNames{feat_idx})(:, 2), '.')
                ylabel(oldNames{feat_idx});
                xlabel('Time (s)');
            end
            
            % Add to feat_idx
            feat_idx = feat_idx + 1;
        end
        
        % Save plot
        saveas(fig, [saveDir, 'Old Feat Figure ', num2str(plot_idx)], ...
            'fig');
        close all
    end
    
else
    % If we have to overlay, then we have two sets of features; otherwise,
    % only one set to plot
    if overlayFeats
        % We assume new and old match up from here on out
        newNames = fieldnames(newFeats);
        oldNames = fieldnames(oldFeats);
        
        newNum = length(newNames); newPlots = ceil(newNum/featPer);
        
        % Take care of plots
        feat_idx = 1;
        for plot_idx = 1:newPlots
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize a tight subplot
            subs = tight_subplot(featPer, 1);
            
            for sub_idx = 1:featPer
                if feat_idx <= newNum
                    % Plot in subplot
                    axes(subs(sub_idx));
                    plot(oldFeats.(oldNames{feat_idx})(:, 1), ...
                        oldFeats.(oldNames{feat_idx})(:, 2), '.')
                    hold on
                    plot(newFeats.(newNames{feat_idx})(:, 1), ...
                        newFeats.(newNames{feat_idx})(:, 2), 'k.')
                    hold off
                    ylabel(newNames{feat_idx});
                    xlabel('Time (s)');
                end
                
                % Add to feat_idx
                feat_idx = feat_idx + 1;
            end
            
            % Save plot
            saveas(fig, [saveDir, 'Overlay Figure ', num2str(plot_idx)], ...
                'fig');
            close all
        end
        
    else
        % We use the single set saved earlier
        singleNum = length(singleNames);
        singlePlots = ceil(singleNum/featPer);
        
        % Take care of single set of features
        feat_idx = 1;
        for plot_idx = 1:singlePlots
            fig = figure(1);
            set(fig, 'Visible', 'on');
            
            % Initialize a tight subplot
            subs = tight_subplot(featPer, 1);
            
            for sub_idx = 1:featPer
                if feat_idx <= singleNum
                    % Plot in subplot
                    axes(subs(sub_idx));
                    plot(singleStruct.(singleNames{feat_idx})(:, 1), ...
                        singleStruct.(singleNames{feat_idx})(:, 2), 'k.')
                    ylabel(singleNames{feat_idx});
                    xlabel('Time (s)');
                end
                
                % Add to feat_idx
                feat_idx = feat_idx + 1;
            end
            
            % Save plot
            saveas(fig, [saveDir, 'Figure ', num2str(plot_idx)], ...
                'fig');
            close all
        end
    end
    
end

end

