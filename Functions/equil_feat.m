function feat_equil = equil_feat(features, method)
% This function takes the features provided to it and equilibrates the
% features according to the method chosen

% ---------- Inputs ------------ %
% features: structure of features stored in Nx2 arrays
% method: string that dictates method used
%         'rest': subtracts away the "rest" value (see code below for
%                 "rest" definition)
%         'mean': subtracts away the mean value

% ------- Outputs -------------- %
% feat_equil: structure of shifted features stored in Nx2 arrays

% If the desired method is to make "rest" the equilibrium
if strcmp(method, 'rest')
    % The amount of time used to calculate "rest" (in seconds)
    restTime = 60;
    
    % Loop through features
    featNames = fieldnames(features);
    for name_idx = 1:numel(featNames)
        % Isolate feature array
        featArray = features.(featNames{name_idx});
        
        % Isolate "rest" data by parsing out the first restTime seconds of
        % data
        restData = featArray(...
            featArray(:, 1) <= featArray(1, 1) + restTime, 2);
        
        % Calculate rest value
        restVal = mean(restData);
        
        % Store and then subtract away
        feat_equil.(featNames{name_idx}) = featArray;
        feat_equil.(featNames{name_idx})(:, 2) = ...
            feat_equil.(featNames{name_idx})(:, 2) - restVal;
    end
    
else
    % Default will be to make the mean the equilibrium
    
    % Loop through features
    featNames = fieldnames(features);
    for name_idx = 1:numel(featNames)
        % Isolate feature array
        featArray = features.(featNames{name_idx});
        
        % Calculate mean value
        meanVal = mean(featArray(:, 2));
        
        % Store and then subtract away
        feat_equil.(featNames{name_idx}) = featArray;
        feat_equil.(featNames{name_idx})(:, 2) = ...
            feat_equil.(featNames{name_idx})(:, 2) - meanVal;
    end
end


end

