function outputStructs = condAvg(inputStructs, varargin)
% This function takes a cell array of structs where responses are separated
% into fields and then takes the repetitions of each response and averages
% to leave just one averaged set of data for each response

% --------- Inputs ---------- %
% inputStructs: cell array of structs, where cell holds one subject's data
%               struct's R fields, where each field contains a TxKxM array,
%               where T is number of timepoints, K is number of features,
%               and M is number of repetitions of the response

% ----- Outputs ---------- %
% outputStructs: cell array of structs, where cell holds one subject's data
%               struct's R fields, where each field contains a TxK array,
%               where T is number of timepoints and K is number of features
% (optional)
% 'firstK_only' - flag that tells the function to only average the first
%                 specified repetitions before storing in outputStructs.
%                 This flag must be followed by a Rx1 or 1xR array of
%                 numbers of repetitions to use for each of the R responses
%                 (i.e., "conditions"). If you want to use all repetitions,
%                 input the number -1. For example, inputting [3, 4, -1]
%                 after 'firstK_only' would average over first 3
%                 repetitions of first response, first 4 repetitions of
%                 second response, and all repetitions of third response

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'firstK_only')
        Kavgs = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('Kavgs', 'var')
    % (arbitarily use the first subject's data)
    % Setting Kavgs to all -1 tells the function to use all repetitions
    Kavgs = -1*ones(length(fieldnames(inputStructs{1})), 1);
end

% Store away fieldnames of the input structs
% (arbitrarily use first subject's data)
Rnames = fieldnames(inputStructs{1});

% Initialize output
outputStructs = cell(length(inputStructs), 1);

% Loop through cells
for sub = 1:length(inputStructs)
    % Loop through responses
    for r = 1:length(Rnames)
        % If we should average over all responses
        if Kavgs(r) == -1
            % Take the mean along the 3rd dimension to average for each
            % response (i.e., "condition") across all repetitions
            % (ignoring NaN because NaN was used as placeholder for missing
            % reptitions' data)
            outputStructs{sub}.(Rnames{r}) = ...
                mean(inputStructs{sub}.(Rnames{r}), ...
                3, 'omitnan');
        else
            % Take the mean along the 3rd dimension for only the
            % repetitions specified
            outputStructs{sub}.(Rnames{r}) = ...
                mean(inputStructs{sub}.(Rnames{r})(:, :, 1:Kavgs(r)), ...
                3, 'omitnan');
        end
    end
end


end