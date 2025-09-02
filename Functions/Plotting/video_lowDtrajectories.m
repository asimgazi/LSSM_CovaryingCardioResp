function video_lowDtrajectories(plotReady, timeVec, dim2plot, saveFile, ...
    varargin)
% This function takes a set of responses separated by protocol
% conditions and participant groups and creates a plot of the specified
% dimensions (2-D or 3-D) such that points of the same timepoint for each
% trajectory are plotted simultaneously frame by frame to create a video
% that can be scrolled through to better understand how the response
% dynamics differ

% ----------- Inputs ------------ %
% plotReady: (# markers) x (# colors) cell array, where each cell
%            contains an array of data to plot corresponding to the color
%            and marker combo specified by the cell array dimensions. These
%            arrays of data will be __xK arrays, where ___ is the variable
%            number of instances and K is the number of dimensions of the
%            response (i.e., # features, principal components, etc.)
% timeVec: Nx1 array storing the first and last time values associated with
%          __ above (i.e., the number of instances). Note that N is equal
%          to the maximum number of timepoints of any of the color-marker
%          combos
% dim2plot: Number of dimensions to plot (2 or 3)
% saveFile: filename to save as (without the .avi at the end)

% (Optional)
% 'plotTitle': flag that is followed by a string that will be used for the
%              title of the figure shown in the video
% 'mrkrNames': flag followed by cell array of labels corresponding to the
%              markers (e.g., groups)
% 'clrNames': flag followed by cell array of labels corresponding to the
%             colors (e.g., conditions)
% 'xlabel': flag followed by x-axis label of the plot
% 'ylabel': flag followed by y-axis label of the plot
% 'zlabel': flag followed by z-axis label of the plot
% 'alpha': flag followed by [alphaMin, alphaMax]
% 'size': flag followed by [sizeMin, sizeMax]
% 'lineWidth': flag followed by [lineWidthMin, lineWidthMax]
% 'velocities': flag that tells the function to create a plot of velocities
%               along each dimension for each of the responses
% 'BSN2023_markers': flag that tells the function to use markers for the
%                    plots according to what I wanted for my BSN 2023
%                    submission where I wanted to use a unique marker for
%                    each repetition, while keeping colors also unique
%                    between conditions and repetitions

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'plotTitle')
        plotTitle = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'mrkrNames')
        mrkrNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'clrNames')
        clrNames = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'xlabel')
        xName = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'ylabel')
        yName = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'zlabel')
        zName = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'alpha')
        alphaMin = varargin{arg + 1}(1);
        alphaMax = varargin{arg + 1}(2);
    elseif strcmp(varargin{arg}, 'size')
        sizeMin = varargin{arg + 1}(1);
        sizeMax = varargin{arg + 1}(2);
    elseif strcmp(varargin{arg}, 'lineWidth')
        lineWidthMin = varargin{arg + 1}(1);
        lineWidthMax = varargin{arg + 1}(2);
    elseif strcmp(varargin{arg}, 'velocities')
        velocityPlots = true;
    elseif strcmp(varargin{arg}, 'BSN2023_markers')
        BSN2023 = true;
    end
end

% Set defaults
if ~exist('plotTitle', 'var'); plotTitle = ''; end
if ~exist('mrkrNames', 'var')
    if size(plotReady, 1) == 1
        % If we only have one "group", no need to even store a name
        mrkrNames = {''};
    else
        mrkrNames = cell(1, size(plotReady, 1));
        for m = 1:length(mrkrNames)
            mrkrNames{m} = ['Group ', num2str(m)];
        end
    end
end
if ~exist('clrNames', 'var')
    clrNames = cell(1, size(plotReady, 2));
    for c = 1:length(clrNames)
        mrkrNames{c} = ['Cond. ', num2str(c)];
    end
end
if ~exist('xName', 'var'); xName = ''; end
if ~exist('yName', 'var'); yName = ''; end
if ~exist('zName', 'var'); zName = ''; end
if ~exist('alphaMin', 'var'); alphaMin = 0.33; end
if ~exist('alphaMax', 'var'); alphaMax = 1.0; end
if ~exist('sizeMin', 'var'); sizeMin = 25; end
if ~exist('sizeMax', 'var'); sizeMax = 75; end
if ~exist('lineWidthMin', 'var'); lineWidthMin = 1.0; end
if ~exist('lineWidthMax', 'var'); lineWidthMax = 3.0; end
if ~exist('velocityPlots', 'var'); velocityPlots = false; end
if ~exist('BSN2023', 'var'); BSN2023 = false; end

% Number of markers needed
numMarkers = size(plotReady, 1);

% Number of colors needed
numColors = size(plotReady, 2);

if numColors <= 3
    % RGB color labels to choose from to differentiate protocol conditions
    % To add more colors for conditions in the future, simply add more rows to
    % this colorLabels array, where each row is a RGB vector
    clrLabels = (1/255)*[221, 85, 255; ...
        255, 85, 85; ...
        95, 211, 188];
elseif numColors > 3 && numColors <= 14
    % If instead we are doing condition separate video plotting, then because
    % the DARPA PTSD dataset has 6 stims, 4 traumas, and 4 neutrals in the day
    % 1 before lunch data, we just hardcode that in here
    
    % Stim. colors
    clrLabels = (1/255)*...
        [141, 95, 211; ...
        212, 42, 255; ...
        255, 85, 221; ...
        68, 33, 120; ...
        222, 135, 205; ...
        85, 0, 212];
    
    % Trauma colors
    clrLabels = [clrLabels; ...
        (1/255)*...
        [255, 42, 127; ...
        255, 85, 85; ...
        211, 95, 95; ...
        211, 95, 141]];
    
    % Neutral colors
    clrLabels = [clrLabels; ...
        (1/255)*...
        [85, 212, 0; ...
        135, 222, 135; ...
        44, 160, 90; ...
        95, 211, 188]];
    
    % If this is for BSN 2023, create associated marker labels, one per
    % condition x repetition combo
    if BSN2023
        mrkrLabels = {'o', 'p', 's', '^', 'd', 'x', ...
            'o', 'p', 's', '^', ...
            'o', 'p', 's', '^'};
    end
else
    % If not enough colors are provided, a random color array will be generated
    rng(5);     % Set random seed so that colors are always the same
    clrLabels = rand(size(plotReady, 2), 3);
end

% Marker labels to choose from to differentiate participant groups
% (to add more markers for groups in the future, simply add more marker
% characters into this 1-D cell array)
if ~BSN2023
    mrkrLabels = {'o', 'p', 's', '^', 'd', 'x', '.', '*', '_', '|', ...
        '+', 'v', '>', '<', 'h'};
end

% If too many groups, throw an error
if numMarkers > numel(mrkrLabels)
    error('Not enough markers available for the number of groups provided');
end

% ------------ Time encoding techniques (besides video) ----------- %
% Encode time in both alpha and size of scatter plot points using
% alpha and size conversion rule below; do the same for line width
timeMin = min(timeVec);
timeMax = max(timeVec);

% Change in transparency and size per time increment
alphaSlope = (alphaMax - alphaMin)/(timeMax - timeMin);
sizeSlope = (sizeMax - sizeMin)/(timeMax - timeMin);
lwSlope = (lineWidthMax- lineWidthMin)/(timeMax - timeMin);

% Store away alpha and size "labels" to encode time; do the same for
% linewidth, except that you need one less element since the lines are
% connecting the markers
alphaLabels = alphaSlope*(timeVec - timeMin) + alphaMin;
sizeLabels = sizeSlope*(timeVec - timeMin) + sizeMin;
lwLabels = lwSlope*(timeVec(2:end) - timeMin) + lineWidthMin;
% ----------------------------------------------------------------- %

% Create figure
fig = figure(1);
set(fig, 'Visible', 'on');  % Make sure figure opens up even in live scripts
set(fig, 'color', 'w') % White background instead of nasty gray
set(fig, 'Position', get(0, 'Screensize'));   % Entire screen
hold on

% Figure out x, y, and z mins and maxes to set plot axes
maxXarray = zeros(numMarkers*numColors, 1);
maxYarray = zeros(numMarkers*numColors, 1);
maxZarray = zeros(numMarkers*numColors, 1);

minXarray = zeros(numMarkers*numColors, 1);
minYarray = zeros(numMarkers*numColors, 1);
minZarray = zeros(numMarkers*numColors, 1);

% Store maxes and mins per response for each dimension
for mrk = 1:numMarkers
    for col = 1:numColors
        % Check for empties (e.g., BSN 2023 causes empty for some conditions)
        if ~isempty(plotReady{mrk, col})
            maxXarray(numColors*(mrk-1) + col) = max(plotReady{mrk, col}(:, 1));
            maxYarray(numColors*(mrk-1) + col) = max(plotReady{mrk, col}(:, 2));
            maxZarray(numColors*(mrk-1) + col) = max(plotReady{mrk, col}(:, 3));
            
            minXarray(numColors*(mrk-1) + col) = min(plotReady{mrk, col}(:, 1));
            minYarray(numColors*(mrk-1) + col) = min(plotReady{mrk, col}(:, 2));
            minZarray(numColors*(mrk-1) + col) = min(plotReady{mrk, col}(:, 3));
        end
    end
end

% Find mins and maxes over all
maxX = (1+sign(max(maxXarray))*.05)*max(maxXarray); 
minX = (1-sign(min(minXarray))*.05)*min(minXarray);
maxY = (1+sign(max(maxYarray))*.05)*max(maxYarray); 
minY = (1-sign(min(minYarray))*.05)*min(minYarray);
maxZ = (1+sign(max(maxZarray))*.05)*max(maxZarray); 
minZ = (1-sign(min(minZarray))*.05)*min(minZarray);

% Branch plotting code based on dimensions to plot (scatter vs. scatter3)
if dim2plot == 2
    % Label axes and title
    % Include variance explained by each of the principal components in the
    % axis labels
    xlabel(xName);
    ylabel(yName);
    title(plotTitle);
    
    % Set x and y limits now so plot doesn't auto rescale in middle of
    % video creation
    xlim([minX, maxX])
    ylim([minY, maxY])
    
    % Create legend using NaN data (so I can label without actually having
    % points show up)
    legendPlots = [];
    legendEntries = cell(numMarkers*numColors, 1);
    for mrk_idx = 1:numMarkers
        for col_idx = 1:numColors
            if ~BSN2023
                L = scatter(nan, nan, sizeMax, clrLabels(col_idx, :), ...
                    mrkrLabels{mrk_idx}, 'filled');
                legendPlots = [legendPlots, L];
                legendEntries{(numColors)*(mrk_idx-1) + col_idx} = ...
                    [mrkrNames{mrk_idx}, '-', clrNames{col_idx}];
            else
                L = scatter(nan, nan, sizeMax, clrLabels(col_idx, :), ...
                    mrkrLabels{col_idx}, 'filled');
                legendPlots = [legendPlots, L];
                legendEntries{(numColors)*(mrk_idx-1) + col_idx} = ...
                    [mrkrNames{mrk_idx}, '-', clrNames{col_idx}];
            end
            
        end
    end
    legend(legendPlots, legendEntries, 'Location', 'bestoutside', ...
        'AutoUpdate', 'off')
    
    % Iterate through timepoints, one at a time, and save away snapshots as
    % video frames
    for t = 1:length(timeVec)
        for mrk_idx = 1:numMarkers
            for col_idx = 1:numColors
                if t < size(plotReady{mrk_idx, col_idx}, 1)
                    % Also plot connecting lines, one edge at a time
                    % Plot lines first so scatter points are on top 
                    % (looks nicer)
                    if t > 1
                        plot([plotReady{mrk_idx, col_idx}(t-1, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 1)], ...
                            [plotReady{mrk_idx, col_idx}(t-1, 2), ...
                            plotReady{mrk_idx, col_idx}(t, 2)], ':', ...
                            'Color', [clrLabels(col_idx, :), alphaLabels(t)], ...
                            'LineWidth', lwLabels(t-1));
                    end
                    
                    % Plot one scatter point at a time to leverage time labels
                    if ~BSN2023
                        s = scatter(plotReady{mrk_idx, col_idx}(t, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 2), sizeLabels(t), ...
                            clrLabels(col_idx, :), mrkrLabels{mrk_idx}, 'filled');
                    else
                        s = scatter(plotReady{mrk_idx, col_idx}(t, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 2), sizeLabels(t), ...
                            clrLabels(col_idx, :), mrkrLabels{col_idx}, 'filled');
                    end
                    
                    % Set the transparency of the scatter point
                    s.MarkerFaceAlpha = alphaLabels(t);
                end
            end
        end
        
        % Store away the video frame for later concatenation into video
        videoFrames(t) = getframe(fig);
        drawnow
    end
    
else
    % Label axes and title
    % Include variance explained by each of the principal components in the
    % axis labels
    xlabel(xName);
    ylabel(yName);
    zlabel(zName);
    title(plotTitle);
    
    % Set x and y limits now so plot doesn't auto rescale in middle of
    % video creation
    xlim([minX, maxX])
    ylim([minY, maxY])
    zlim([minZ, maxZ])
    
    % Create legend using NaN data (so I can label without actually having
    % points show up)
    % Create legend using NaN data (so I can label without actually having
    % points show up)
    legendPlots = [];
    legendEntries = cell(numMarkers*numColors, 1);
    for mrk_idx = 1:numMarkers
        for col_idx = 1:numColors
            if ~BSN2023
                L = scatter3(nan, nan, nan, sizeMax, clrLabels(col_idx, :), ...
                    mrkrLabels{mrk_idx}, 'filled');
                legendPlots = [legendPlots, L];
                legendEntries{(numColors)*(mrk_idx-1) + col_idx} = ...
                    [mrkrNames{mrk_idx}, '-', clrNames{col_idx}];
            else
                L = scatter3(nan, nan, nan, sizeMax, clrLabels(col_idx, :), ...
                    mrkrLabels{col_idx}, 'filled');
                legendPlots = [legendPlots, L];
                legendEntries{(numColors)*(mrk_idx-1) + col_idx} = ...
                    [mrkrNames{mrk_idx}, '-', clrNames{col_idx}];
            end
        end
    end
    legend(legendPlots, legendEntries, 'Location', 'bestoutside', ...
        'AutoUpdate', 'off')
    
    % Iterate through timepoints, one at a time, and save away snapshots as
    % video frames
    for t = 1:length(timeVec)
        for mrk_idx = 1:numMarkers
            for col_idx = 1:numColors
                if t < size(plotReady{mrk_idx, col_idx}, 1)
                    % Also plot connecting lines, one edge at a time
                    % I plot line first, so markers are on top (look nicer)
                    if t > 1
                        plot3([plotReady{mrk_idx, col_idx}(t-1, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 1)], ...
                            [plotReady{mrk_idx, col_idx}(t-1, 2), ...
                            plotReady{mrk_idx, col_idx}(t, 2)], ...
                            [plotReady{mrk_idx, col_idx}(t-1, 3), ...
                            plotReady{mrk_idx, col_idx}(t, 3)], ':', ...
                            'Color', [clrLabels(col_idx, :), alphaLabels(t)], ...
                            'LineWidth', lwLabels(t-1));
                    end
                    
                    % Plot one scatter point at a time to leverage time labels
                    if ~BSN2023
                        s = scatter3(plotReady{mrk_idx, col_idx}(t, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 2), ...
                            plotReady{mrk_idx, col_idx}(t, 3), ...
                            sizeLabels(t), ...
                            clrLabels(col_idx, :), mrkrLabels{mrk_idx}, 'filled');
                    else
                        s = scatter3(plotReady{mrk_idx, col_idx}(t, 1), ...
                            plotReady{mrk_idx, col_idx}(t, 2), ...
                            plotReady{mrk_idx, col_idx}(t, 3), ...
                            sizeLabels(t), ...
                            clrLabels(col_idx, :), mrkrLabels{col_idx}, 'filled');
                    end
                    
                    % Set the transparency of the scatter point
                    s.MarkerFaceAlpha = alphaLabels(t);
                end
            end
        end
        
        % Store away the video frame for later concatenation into video
        videoFrames(t) = getframe(fig);
        drawnow
    end
end

hold off

% Store one final frame to capture the legend
videoFrames(end+1) = getframe(fig);
drawnow

% Save final frame as its own separate static figure
saveas(fig, saveFile, 'fig');

% Create a video object
writerObj = VideoWriter([saveFile, '.avi']);
writerObj.FrameRate = 7;   % Frames shown per second in video

% Open the video writer
open(writerObj);

% Check the frame dimensions of the first frame so I can resize if ever
% needed to avoid errors (sometimes, randomly one of the frames will be off
% by like one pixel in either height or width and cause code to error)
frameHeight = size(videoFrames(1).cdata, 1);
frameWidth = size(videoFrames(1).cdata, 2);

% Write the frames to the video
for f = 1:length(videoFrames)
    frame = videoFrames(f);
    
    % Resize if necessary
    if size(frame.cdata, 1) ~= frameHeight || ...
            size(frame.cdata, 2) ~= frameWidth
        frame = imresize(frame.cdata, [frameHeight, frameWidth]);
    end
    
    writeVideo(writerObj, frame);
end

% Close the writer object
close(writerObj);

% If the user wanted to create velocity plots, do so now
if velocityPlots
    % Will store velocities
    velocs = cell(numMarkers, numColors);
    
    % Just to make sure code below works, we need to make sure timeVec is a
    % column vector
    if size(timeVec, 1) < size(timeVec, 2)
        timeVec = timeVec';
    end
    
    % Create figure
    fig = figure();
    set(fig, 'Visible', 'on');  % Make sure figure opens up even in live scripts
    
    % Loop through dimensions
    for dim = 1:dim2plot
        subplot(dim2plot + 1, 1, dim);
        hold on
        % Loop through markers and colors
        for mrk_idx = 1:numMarkers
            for col_idx = 1:numColors
                % Compute velocity for this response along this dimension
                velocs{mrk_idx, col_idx}(:, dim) = ...
                    diff(plotReady{mrk_idx, col_idx}(:, dim))./...
                    diff(timeVec(1:size(plotReady{mrk_idx, col_idx}, 1)));
                
                % Create a simple plot of velocity vs. time
                plot(timeVec(1:(size(plotReady{mrk_idx, col_idx}, 1)-1)), ...
                    velocs{mrk_idx, col_idx}(:, dim), 'Marker', mrkrLabels{mrk_idx}, ...
                    'Color', clrLabels(col_idx, :), ...
                    'MarkerEdgeColor', clrLabels(col_idx, :), ...
                    'MarkerFaceColor', clrLabels(col_idx, :));
            end
        end
        % Create a legend using the legend entries already created
        legend(legendEntries);
        title(['Velocity Along Dim. ', num2str(dim)])
        xlabel('Time (s)')
        ylabel('Velocity (\sigma/s)')
        hold off
    end
    
    % Take care of speed
    speed = cell(numMarkers, numColors);
    subplot(dim2plot + 1, 1, dim2plot + 1);
    hold on
    for mrk_idx = 1:numMarkers
        for col_idx = 1:numColors
            % Compute speed
            for t = 1:size(velocs{mrk_idx, col_idx}, 1)
                speed{mrk_idx, col_idx}(t) = norm(velocs{mrk_idx, col_idx}(t, :));
            end
            
            % Plot speed
            plot(timeVec(1:(size(plotReady{mrk_idx, col_idx}, 1)-1)), ...
                speed{mrk_idx, col_idx}, 'Marker', mrkrLabels{mrk_idx}, ...
                'Color', clrLabels(col_idx, :), ...
                'MarkerEdgeColor', clrLabels(col_idx, :), ...
                'MarkerFaceColor', clrLabels(col_idx, :));
        end
    end
    % Create a legend using the legend entries already created
    legend(legendEntries);
    title('Speed')
    xlabel('Time (s)')
    ylabel('Speed (\sigma/s)')
    hold off
    
    % Title the whole subplot grid
    sgtitle('Velocities and Speed')
    
    % Save the figure
    saveas(fig, [saveFile, '_velocities'], 'fig');
end

end

