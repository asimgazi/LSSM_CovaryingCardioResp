function resampled_featStruct = resamp_feat(featStruct,resamp_Fs, ...
    resamp_method)
% This function takes a structure of features and resamples each of the
% features to the desired frequency using the desired resampling method

% ---------- Inputs --------------- %
% featStruct: struct containing a feature per field, where each field is a
%             __x2 array, with the first column being time and the second 
%             column holding the feature values
% resamp_Fs: frequency to resample features to
% resamp_method: resampling method to use in the clean_resample.m function

% -------- Outputs ------------ %
% resampled_featStruct: struct containing a feature per field, where each 
%                       field now contains a ___x2 array with the first
%                       column holding time where intersample intervals are
%                       now dictated by resamp_Fs and the second column
%                       holds the corresponding resampled values

% Get the feature names to loop through
featNames = fieldnames(featStruct);

% Loop through feature names
for name_idx = 1:numel(featNames)
    % An issue we run into when resampling and trying to overwrite the
    % same data structure with the output is that we have
    % to do so each column separately. The resampled
    % output, however, is of a different number of rows.
    % So we will need to use a dummy vector
    dummyArr = [];
    
    % Resample and store in dummy array
    [dummyArr(:, 2), dummyArr(:, 1)] = ...
        clean_resample(featStruct.(featNames{name_idx})(:, 2), ...
        featStruct.(featNames{name_idx})(:, 1), resamp_Fs, resamp_method);
    
    % Store away in final result
    resampled_featStruct.(featNames{name_idx}) = dummyArr;
end

end

