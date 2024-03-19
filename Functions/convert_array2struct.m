function outputStructs = convert_array2struct(inputArrays, condNames)
% This function takes a cell array of arrays formatted such that you have
% multiple responses of multiple features for a fixed amount of time and
% generalizes into a struct that separates the responses into their own
% fields such that they can each be named and can differ in time length

% --------- Inputs ---------- %
% inputArrays: cell array of arrays, where cell holds one subject's data
%             Array is TxMxR, where T is number of time points, M is number
%             of features, and R is number of responses (responses will
%             eventually be stored separately)
% condNames: cell array of length R that stores a string name for each
%            response to be used as fieldnames of the output struct. NOTE:
%            the ordering of these names needs to match the ordering of
%            data along the 3rd dimension (i.e., R dimension) of arrays

% ----- Outputs ---------- %
% outputStructs: cell array of structs, where each cell holds one subject's
%               data. The struct has R fields, where each field stores a
%               TxM array

% Quick dummy proofing check
for sub = 1:length(inputArrays)
    if length(condNames) ~= size(inputArrays{sub}, 3)
        error(['Things are not lining up for cell number ', num2str(sub)])
    end
end

% Initialize output based on length of input
outputStructs = cell(length(inputArrays), 1);

% Loop through length of cell array and store away data
for sub = 1:length(inputArrays)
    % Initialize this subject's struct
    outputStructs{sub} = struct();
    
    % Loop through responses
    for resp = 1:length(condNames)
        % Store away as desired
        outputStructs{sub}.(condNames{resp}) = inputArrays{sub}(:, :, resp);
    end
end


end