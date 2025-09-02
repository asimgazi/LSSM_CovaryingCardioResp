function fitModel = ...
    fit_LSSM(timeVec, featArray, modelOrder, ...
    inputArray, inputOption, varargin)
% This function takes time, features, inputs, and hyperparameters and fits
% a linear state-space model (LSSM)

% --------- Inputs ----------- %
% timeVec: Nx1 array of timepoints
% featArray: NxM array of M features for N timepoints
% modelOrder: model order to use
% inputArray: NxP array of P inputs for N timepoints
% inputOption: string that dictates what inputs were included and what
%              delay options are to come in varargin
%              'v': VNS only
%              'vt': VNS and trauma recall
%              'vn': VNS and neutral
%              'vtn': VNS, trauma recall, and neutral
% varargin: cell array that will contain the corresponding input delays to
%           when fitting the model

% -------- Outputs ---------- %
% fitModel: the LSSM model fit using the parameters input to the function


% Create an iddata object using the features and inputs
% (data are uniformly sampled, so I can just take intersample interval for
% first two elements to get Ts)
sampleTime = timeVec(2) - timeVec(1);
fitData = iddata(featArray, inputArray, sampleTime);

% Enforce stability for ssest estimation
% Stability ensures that states return to equilibrium once inputs subside
ssOptions = ssestOptions('EnforceStability', true);

% Proceed according to inputOption
if strcmp(inputOption, 'v')
    vnsDelay = varargin{1};
    
    % Fit the model using ssest()
    % N4SID estimation followed by prediction error minimization
    fitModel = ssest(fitData, modelOrder, ssOptions, ...
        'InputDelay', vnsDelay, 'Form', 'modal', 'Ts', sampleTime);
    
elseif strcmp(inputOption, 'vt')
    vnsDelay = varargin{1};
    traumDelay = varargin{2};
    
    % Fit the model using ssest()
    % N4SID estimation followed by prediction error minimization
    fitModel = ssest(fitData, modelOrder, ssOptions, ...
        'InputDelay', [vnsDelay; traumDelay], ...
        'Form', 'modal', 'Ts', sampleTime);    
    
elseif strcmp(inputOption, 'vn')
    vnsDelay = varargin{1};
    neutrDelay = varargin{2};
    
    % Fit the model using ssest()
    % N4SID estimation followed by prediction error minimization
    fitModel = ssest(fitData, modelOrder, ssOptions, ...
        'InputDelay', [vnsDelay; neutrDelay], ...
        'Form', 'modal', 'Ts', sampleTime);
    
else
    vnsDelay = varargin{1};
    traumDelay = varargin{2};
    neutrDelay = varargin{3};
    
    % Fit the model using ssest()
    % N4SID estimation followed by prediction error minimization
    fitModel = ssest(fitData, modelOrder, ssOptions, ...
        'InputDelay', [vnsDelay; traumDelay; neutrDelay], ...
        'Form', 'modal', 'Ts', sampleTime);
end

