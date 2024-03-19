function plotArrays(featPer, saveDir, varargin)
% This function takes structs and plots them according to the options
% provided

% ------ Input Instructions ------- %
% featPer: features per plot
% saveDir: the directory to save the plots to
% 1. One of the inputs needs to be a string that is either 'overlay' or
%    'separate'. This string will dictate whether the columns of arrays
%    are overlaid on top of each other or kept in separate figures. By
%    default, each feature within will be on a separate subplot,
%    but this option is more for a new set of features and an old set of
%    features, for example. This can be used to understand the effects of a
%    particular operation on the features
% 2. One of the inputs needs to be a time vector that is the same length as
%    arrays input to the function. This time vector needs to be input after
%    the string 'time'. If the new and old arrays have different time
%    vectors, input the time vector associated with new data after
%    'timeNew' and the time vector associated with old data after 'timeOld'
% 3. If you want to input an array of features or inputs, you need to 
%    either specify whether the set of features is 'new' or 'old'. 
%    This function, as of 1/9/2022, allows for two arrays max, where these
%    features can either be plotted separately or overlaid. To do this, you
%    first input the string 'new' and then you input the features next. You
%    also need to input a cell array of names afterwards

% -------- Use examples -------- %
% plotArrays(6, ['Test Figures', filesep], 'overlay', ...
% 'new', newFeats, newFeatNames, ...
% 'old', oldFeats, oldFeatNames, 'time', timeVec);
% plotArrays(12, ['Wow Figures', filesep], 'separate', ...
% 'new', feats, featNames, 'time', timeVec);


% We first need to parse the varargin looking for whether to overlay
% features if provided with more than one struct
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'overlay')
        overlay = true;
    elseif strcmp(varargin{i}, 'new')
        newArray = varargin{i+1};
        newNames = varargin{i+2};
    elseif strcmp(varargin{i}, 'old')
        oldArray = varargin{i+1};
        oldNames = varargin{i+2};
    elseif strcmp(varargin{i}, 'time')
        timeVec = varargin{i+1};
    elseif strcmp(varargin{i}, 'timeNew')
        timeNew = varargin{i+1};
    elseif strcmp(varargin{i}, 'timeOld')
        timeOld = varargin{i+1};
    end
end

% Set defaults and booleans to know how to proceed later
if ~exist('overlay', 'var')
    overlay = false;
end
if overlay && exist('newArray', 'var') && exist('oldArray', 'var')
    twoPlots = false;
elseif exist('newArray', 'var') && exist('oldArray', 'var')
    twoPlots = true;
else
    twoPlots = false;
    
    % Assign the single array so I can easily work with it
    % later without checking for which one exists
    if exist('newArray', 'var')
        singleTime = timeVec;
        singleArray = newArray;
        singleNames = newNames;
    else
        singleTime = timeVec;
        singleArray = oldArray;
        singleNames = oldNames;
    end
end

% Figure out if time vector is shared
if ~exist('timeNew', 'var')
    timeNew = timeVec;
    timeOld = timeVec;
end


% We now know how to proceed! :)


% Split procedures based on two plots or one
if twoPlots
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
                plot(timeNew, newArray(:, feat_idx), 'k.')
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
                plot(timeOld, oldArray(:, feat_idx), '.')
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
    if overlay
        % Figure out how many we're dealing with
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
                    plot(timeOld, oldArray(:, feat_idx), '.')
                    hold on
                    plot(timeNew, newArray(:, feat_idx), 'k.')
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
                    plot(singleTime, singleArray(:, feat_idx), 'k.')
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

