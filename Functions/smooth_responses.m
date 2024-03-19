function smoothResps = smooth_responses(inputResps, smoothTime, Fs, varargin)
% This function takes a "feature array" and smooths each feature according
% to some time window to smooth over

% -------- Inputs ----------- %
% inputResps: NxK array, where N represents number of samples of data
%             available and K represents the dimensionality of responses 
%             (i.e., "features")
% smoothTime: the time window to smooth over
% Fs: sampling frequency of the data
% (Optional)
% 'method': If you input this flag, then the next input string should be
%           the method input to the smooth() MATLAB function

% ----------- Outputs ---------- %
% smoothResps: NxK array storing the smoothed version of the input

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'method')
        smoothMethod = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('smoothMethod', 'var'); smoothMethod = 'moving'; end

% Compute smoothing "span" (i.e., number of samples to smooth over) using
% sampling frequency and smoothing time window
smoothSpan = round(smoothTime*Fs);

% Initialize smoothed "features" array
smoothResps = zeros(size(inputResps));

% Loop through "features" and smooth each one
for feat_idx = 1:size(inputResps, 2)
    smoothResps(:, feat_idx) = smooth(inputResps(:, feat_idx), ...
        smoothSpan, smoothMethod);
end


end

