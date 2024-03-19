function naivePreds = naivePred(trueData, varargin)
% This function takes a vector of true data and simply returns a naive
% predictor's corresponding predictions

% ----- Inputs ------- %
% trueData: NxK array of true datapoints (e.g., test data), where
%           N is number of instances and K is number of features/outputs
% (optional)
% 'steps' - flag that is followed by number of steps ahead the naive
% predictor must predict (default if not specified is steps = 1)
% 'initialEst' - flag that is followed by an initial prediction to make
% (default if not specified is y = 0) --> Note that the initialEst
% vector must be either of dimension steps x K or 1 x K. If 1 x K, then the
% row vector will be copied steps times to fill all initial output
% estimates

% ----- Outputs -------- %
% naivePreds: NxK vector of naive predictions

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'steps')
        steps = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'initialEst')
        initialEst = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('steps', 'var')
    steps = 1;
end
if ~exist('initialEst', 'var')
    initialEst = zeros(steps, size(trueData, 2));
end

% If only one initial estimate provided but steps > 1, then we'll just copy
% that estimate steps times
if size(initialEst, 1) == 1 && steps > 1
    initialEst(1:steps, :) = repmat(initialEst, steps, 1);
end

% Initialize naivePreds array
naivePreds = zeros(size(trueData));

% Loop through columns
for col = 1:size(trueData, 2)
    % Store initial condition of naive prediction
    naivePreds(1:steps, col) = initialEst(:, col);
    
    % The naive predictor simply predicts the next time point will be the
    % time point K before, i.e., naivePreds(K + 1) = trueData(1),
    % naivePreds(K + 2) = trueData(2), ..., 
    % naivePreds(K + N) = trueData(N)
    for row = (steps+1):size(trueData, 1)
        naivePreds(row, col) = trueData(row-steps, col);
    end
end

end

