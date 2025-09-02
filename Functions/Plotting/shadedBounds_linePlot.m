function shadedBounds_linePlot(center_TS, upperBound_TS, lowerBound_TS, ...
    timeVec, color, varargin)
% This function takes a center time series, an upper bound time series, and
% a lower bound time series and creates line plots with shaded bounds to
% depict upper and lower bounds (e.g., mean +- SEM plots)
% Code adapted from some example code I found online

% -------- Inputs ------------ %
% center_TS: 1xT or Tx1 vector of time points of the center time series
%            that will be the solid line
% upperBound_TS: 1xT or Tx1 vector of time points of the upper bound time
%                series that will serve as upper bound of shaded region
% lowerBound_TS: 1xT or Tx1 vector of time points of the lower bound time
%                series that will serve as lower bound of shaded region
% timeVec: 1XT or Tx1 time vector to plot against
% color: either RGB triplet or named color to color the plot using
% (Optional)
% you can input a marker character to use for line plot (e.g., 'o')

% Make everything row vectors so we do not have problems when working with
% fliplr() function to create polygon to fill in with a patch

if ~isempty(varargin)
    % Store the marker character
    markerChar = varargin{1};
else
    markerChar = '.';
end

if size(center_TS, 2) == 1
    center_TS = center_TS';
end
if size(upperBound_TS, 2) == 1
    upperBound_TS = upperBound_TS';
end
if size(lowerBound_TS, 2) == 1
    lowerBound_TS = lowerBound_TS';
end
if size(timeVec, 2) == 1
    timeVec = timeVec';
end

% Create a time vector that goes forward and backward for filling between
t2 = [timeVec, fliplr(timeVec)];

% Create corresponding upper and lower bounds to fill in between
inBetween = [upperBound_TS, fliplr(lowerBound_TS)];

hold on
% Fill in between and set shading (also don't let it show up in the legend)
patch(t2, inBetween, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

% Plot the center time series on top and include in legend
plot(timeVec, center_TS, ['-', markerChar], 'Color', color, ...
    'LineWidth', 0.7, 'MarkerFaceColor', color, 'MarkerSize', 3);
hold off
end

