function alignedFeat = align_feat(resampFeat, resamp_Fs, interpMethod)
% This function takes features resampled to same Fs and aligns them 
% according to the last start time and first end time of all the features

% -------- Inputs ------------- %
% resampFeat: structure of resampled features, where each field stores a
%             ____x2 array of time (first column) and data (second column) 
%             for one feature
% resamp_Fs: frequency that resampFeat were resampled to
% interpMethod: interpolation method used in interp1 function


% -------- Outputs ----------- %
% alignedFeat: structure of resampled and aligned features, where each
%              field stores a Nx2 array of time (first column) and data 
%              (second column) for a feature (notice that N is now fixed
%              because the lengths of all features are now the same)

% Feature names to loop through
featNames = fieldnames(resampFeat);

% Initialize start and end times for all of the features
startTimes = zeros(numel(featNames), 1);
endTimes = zeros(numel(featNames), 1);

% Store start and end times for all of the features
for name_idx = 1:numel(featNames)
    startTimes(name_idx) = resampFeat.(featNames{name_idx})(1, 1);
    endTimes(name_idx) = resampFeat.(featNames{name_idx})(end, 1);
end

% Find max start and min end times
maxStart = max(startTimes);
minEnd = min(endTimes);

% Create new time vector to query according to
newTime = (maxStart:(1/resamp_Fs):minEnd)';

% Use interp1 to resize, leveraging a dummy array for storing purposes
for name_idx = 1:numel(featNames)
    dummyArr = zeros(length(newTime), 2);
    dummyArr(:, 1) = newTime;
    dummyArr(:, 2) = interp1(resampFeat.(featNames{name_idx})(:, 1), ...
        resampFeat.(featNames{name_idx})(:, 2), newTime, interpMethod);
    
    % Store away final result
    alignedFeat.(featNames{name_idx}) = dummyArr;
end

end

