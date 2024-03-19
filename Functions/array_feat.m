function [featArray, timeVec, featNames] = array_feat(featStruct)
% This function takes features stored in struct form and returns instead an
% array of features, wehre each feature is a column in the order
% corresponding to featNames. A timeVec is also returned to keep track of
% time

% -------- Inputs ----------- %
% featStruct: structure of features, where each field stores an Nx2 array

% ------- Outputs ----------- %
% featArray: N x (# of features) array of feature values
% timeVec: N x 1 array storing corresponding time values
% featNames: cell array with (# of features) cells, where each cell stores
%            a feature name, in the same order as the columns of featArray

% Get the names of the features
featNames = fieldnames(featStruct);

% Initialize feature array (all features are same length)
featArray = zeros(size(featStruct.(featNames{1}), 1), numel(featNames));

% Loop through features and store away feature values in each column
for name_idx = 1:numel(featNames)
    featArray(:, name_idx) = featStruct.(featNames{name_idx})(:, 2);
end

% Store away time
timeVec = featStruct.(featNames{1})(:, 1);

end

