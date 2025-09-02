function [feat_desir, featNames_desir] = combo_desirFeat(varargin)

% This function takes in feature structures and the names of the desired
% features and combines them into a single struct storing all the desired
% features

% ---------- Inputs ------------ %
% varargin: Inputs should come in pairs, where first of each pair is a 1x2
%           cell array, where each cell corresponds to before/after lunch
%           struct of features. Each struct field should be a __ x 2 array,
%           first column time, second column data; UNLESS, the feature is a
%           window-based feature. In that case, ___ x 3 array, where first
%           column is window start time, second column is window end time,
%           and third column is data
%           Second in each pair should be the cell array of corresponding
%           names desired to extract from the struct of features
% Example function call:
%  [feat_desir, featNames_desir] = combo_desirFeat(beatFeatures_struct,
%  beatFeatures_names, respFeatures_struct, respFeatures_names)

% -------- Outputs ----------- %
% feat_desir:   Struct storing all desired features, where each field is
%               _x2 - first column is time, second column is data
% featNames_desir: Feature names in a cell array in order


% Initialize final result
feat_desir = struct();

% Loop through the varargin in pairs - WE WANT BEFORE LUNCH ONLY
for in_idx = 1:2:numel(varargin)
    % Loop through names
    for name_idx = 1:numel(varargin{in_idx + 1})
        % If the array has 3 columns, then we have a window-based feature
        if size(varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx}), ...
                2) > 2
            % Initialize 2 column array for final result
            feat_desir.(varargin{in_idx + 1}{name_idx}) = zeros(size(...
                varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx}), ...
                1), 2);
            
            % 1st column = time
            % Our approach, as of 12/26/2021, is to use the window midpoint
            feat_desir.(varargin{in_idx + 1}{name_idx})(:, 1) = ...
                0.5*(varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx}...
                )(:, 1) + ...
                varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx})...
                (:, 2));
            
            % 2nd column = data
            feat_desir.(varargin{in_idx + 1}{name_idx})(:, 2) = ...
                varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx})(:, 3);
        else
            % Otherwise, we can simply assign and move on
           feat_desir.(varargin{in_idx + 1}{name_idx}) = ...
                varargin{in_idx}{1}.(varargin{in_idx + 1}{name_idx}); 
        end
    end
end

% Store away all of the field names in order
featNames_desir = fieldnames(feat_desir);

end

