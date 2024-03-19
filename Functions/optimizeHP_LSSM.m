function [mOrder_star, varargout] = ...
    optimizeHP_LSSM(verbose, timeVec, featArray, modelOrders, ...
    inputArray, inputOption, varargin)
% This function takes time, features, and inputs, along with
% hyperparameters to optimize over, and returns the optimal hyperparameters
% according to some model fit metric

% (I later generalized function to deal with generic "loss" rather than
% specifically AICc, so comments may not always agree in that regard)

% --------- Inputs ----------- %
% verbose: boolean to decide whether to output model estimation status
% timeVec: Nx1 array of timepoints
% featArray: NxM array of M features for N timepoints
% modelOrders: 1-D array of latent state dimensionality to optimize over
% inputArray: NxP array of P inputs for N timepoints
% inputOption: string that dictates what inputs were included and what
%              delay options are to come in varargin
%              'v': VNS only
%              'vt': VNS and trauma recall
%              'vn': VNS and neutral
%              'vtn': VNS, trauma recall, and neutral
% varargin: cell array that will contain the corresponding input delays to
%           optimize over based on inputOptions (e.g., if inputOption is
%           set to 'vt', then varargin{1} will be the set of VNS delays to
%           optimize over, varargin{2} will be the set of trauma delays
% (Optional)
% 'max_mOrder': flag that tells this function to just use the maximum model
%               order allowable (i.e., that does not result in an improper
%               fit). NOTE: this flag must be provided AFTER providing the
%               delays specified above (see varargin)
% 'leak_test': flag that tells this function test data will be
%              provided to use to optimize hyperparameters (i.e.,
%              purposeful data leakage where model will still be fit using
%              training data, but testing data will be used for evaluation
%              to optimize model order). If this flag is provided, then the
%              next three inputs need to be timeTest, featTest, and
%              inputTest (i.e., time, features, and inputs for test data in
%              same format as train data, except of course with different
%              number of rows)
% 'validation': flag that tells this function that a validation data
%               approach will be used to optimize hyperparameters (i.e., a
%               section of the training data will be allocated for
%               validation). If this flag is provided, then the next input
%               must be a fractional amount of the final portion of the
%               training data to be allocated for validation. Percentages
%               can also be provided (i.e., 20 or 0.2 produce same result)

% -------- Outputs ---------- %
% mOrder_star: optimal model order (dimension of latent state), either
%              based on optimization routine or flag specification(s) 
%              (see optional inputs)
% varargout: cell array that will contain the optimal input delays
%            corresponding to varargin

% Check for flags in varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'max_mOrder')
        % Set this boolean to true to maximize model order
        max_mOrder = true;
    elseif strcmp(varargin{arg}, 'leak_test')
        % Set this boolean to true for purposeful data leakage
        dataLeak = true;
        
        % Store test data
        timeTest = varargin{arg + 1};
        featTest = varargin{arg + 2};
        inputTest = varargin{arg + 3};
    elseif strcmp(varargin{arg}, 'validation')
        % Set this boolean to true for validation data for optimization
        valApproach = true;
        
        % Handle whether user gives percentage or decimal
        if varargin{arg + 1} > 1
            prctVal = varargin{arg + 1}/100;
        else
            prctVal = varargin{arg + 1};
        end
        
        % Figure out how many validation data points to set aside
        numVal = round(prctVal*length(timeVec));
        
        % Separate validation data from pure training
        timeVal = timeVec((end - numVal + 1):end);
        featVal = featArray((end - numVal + 1):end, :);
        inputVal = inputArray((end - numVal + 1):end, :);
        
        % Reassign so reuse of code is easy
        timeVec = timeVec(1:(end-numVal));
        featArray = featArray(1:(end-numVal), :);
        inputArray = inputArray(1:(end-numVal), :);
    end
end

% Set defaults
if ~exist('max_mOrder', 'var'); max_mOrder = false; end
if ~exist('dataLeak', 'var'); dataLeak = false; end
if ~exist('valData', 'var'); valApproach = false; end


% We suppress these expected warnings if verbosity is not desired
% (we already ignore these models when optimizing HPs, so we might as well
% ignore the warnings)
if ~verbose
    % Specify the warning IDs to turn off
    wrnID = {'Ident:estimation:NparGTNsamp', ...
        'Ident:estimation:illConditionedCovar2'};
    
    % Turn these warnings off
    for id = 1:numel(wrnID)
        warning('off', wrnID{id})
    end
end

% Create an iddata object using the features and inputs
% (data are uniformly sampled, so I can just take intersample interval for
% first two elements to get Ts)
sampleTime = timeVec(2) - timeVec(1);
fitData = iddata(featArray, inputArray, sampleTime);

% If we are going to use test data to optimize, create iddata object with
% test data
if dataLeak
    testData = iddata(featTest, inputTest, sampleTime);
elseif valApproach
    % If we are going to use validation data to optimize, create iddata
    % object with validation data
    valData = iddata(featVal, inputVal, sampleTime);
end

% Enforce stability for later ssest estimation
% Stability ensures that states return to equilibrium once inputs subside
ssOptions = ssestOptions('EnforceStability', true);

% Proceed according to inputOption
if strcmp(inputOption, 'v')
    varargout = cell(1, 1);
    vnsDelays = varargin{1};
    
    % Check if we're supposed to just use the maximum available
    if max_mOrder
        % Initialize array to store AICc, but only allow one option for
        % model order
        model_loss = zeros(1, length(vnsDelays));
        
        % Since we need to use the maximum available model order, we should
        % loop backwards until we fit the model error free
        m_idx = length(modelOrders);    % Initialize at max
        fitSuccess = false;             % Initialize boolean to false
        while ~fitSuccess
            % Start with true for looping through delays to test out
            loopDelays = true; stop_m_idx = m_idx;
            for v_idx = 1:length(vnsDelays)
                if loopDelays
                    % Since input delay should have no effect on whether the
                    % model was successfully fit (should have no effect on
                    % number of parameters, which means no consideration of
                    % number of datapoints vs. number of parameters), once the
                    % model was successfully fit once, then we can assume it
                    % will be for all possible delays, so we set fitSuccess to
                    % true, never check again, and start storing away AICc
                    
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', vDelay, 'Form', 'modal', 'Ts', sampleTime);
                    
                    % If the number of free parameters is less than the number
                    % of datapoints available, then I skip this model order
                    if isempty(modelFit.Report.Fit.AICc)
                        % Set looping boolean to false so we skip the rest of
                        % the delays
                        loopDelays = false;
                    else
                        % If we were successful, then it's time to just stick
                        % with this model order and loop through delays
                        fitSuccess = true;
                    end
                    
                    % If the covariance matrix estimate is unreliable, then I
                    % set the AICc value to infinity so that model does not get
                    % selected (only necessary if loopDelays is still true)
                    if loopDelays && isempty(modelFit.covar)
                        model_loss(1, v_idx) = inf;
                    elseif loopDelays
                        % Store AICc
                        model_loss(1, v_idx) = ...
                            modelFit.Report.Fit.AICc;
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Max model order mode: Done with model order ', ...
                            num2str(mOrder), ' and VNS delay ', ...
                            num2str(vDelay), ' at ', datestr(timeNow)])
                    end
                end
            end
            % Negative increment
            m_idx = m_idx - 1;
        end
        
        % Check if we need to purposefully leak data for HP optimization
    elseif dataLeak
        % Initialize array to store prediction errors on test data for each
        % model
        model_loss = zeros(length(modelOrders), length(vnsDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                % Model order to use in this iteration
                mOrder = modelOrders(m_idx);
                
                % VNS delay to use in this iteration
                vDelay =  vnsDelays(v_idx);
                
                % Estimate models using ssest()
                % N4SID estimation followed by prediction error minimization
                modelFit = ssest(fitData, mOrder, ssOptions, ...
                    'InputDelay', vDelay, 'Form', 'modal', 'Ts', sampleTime);
                
                % If the covariance matrix estimate is unreliable or
                % the number of free parameters is less than the number
                % of datapoints available, then I set the loss value to
                % infinity so that model does not get selected
                if isempty(modelFit.covar) || ...
                        isempty(modelFit.Report.Fit.AICc)
                    model_loss(m_idx, v_idx) = inf;
                else
                    % Evaluate model prediction on test data
                    R2_test = eval_modelPred(modelFit, testData, ...
                        'K_steps', 1, 'Rsquared', ...
                        'IC_training', inputArray, featArray);
                    
                    % Store loss on test data (loss will just be negative
                    % of the mean R^2 that we use when evaluating later)
                    model_loss(m_idx, v_idx) = -1*mean(R2_test, 2);
                end
                
                % Display model estimation status for user
                if verbose
                    timeNow = datetime('now','TimeZone','local',...
                        'Format','d-MMM-y HH:mm:ss Z');
                    disp(['Done with model order ', num2str(mOrder), ...
                        ' and VNS delay ', num2str(vDelay), ...
                        ' at ', datestr(timeNow)])
                end
            end
        end
        
    % Check if we need to use validation data for HP optimization
    elseif valApproach
        % Initialize array to store prediction errors on val. data for each
        % model
        model_loss = zeros(length(modelOrders), length(vnsDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                % Model order to use in this iteration
                mOrder = modelOrders(m_idx);
                
                % VNS delay to use in this iteration
                vDelay =  vnsDelays(v_idx);
                
                % Estimate models using ssest()
                % N4SID estimation followed by prediction error minimization
                modelFit = ssest(fitData, mOrder, ssOptions, ...
                    'InputDelay', vDelay, 'Form', 'modal', 'Ts', sampleTime);
                
                % If the covariance matrix estimate is unreliable or
                % the number of free parameters is less than the number
                % of datapoints available, then I set the loss value to
                % infinity so that model does not get selected
                if isempty(modelFit.covar) || ...
                        isempty(modelFit.Report.Fit.AICc)
                    model_loss(m_idx, v_idx) = inf;
                else
                    % Evaluate model prediction on validation data
                    R2_val = eval_modelPred(modelFit, valData, ...
                        'K_steps', 1, 'Rsquared', ...
                        'IC_training', inputArray, featArray);
                    
                    % Store loss on test data (loss will just be negative
                    % of the mean R^2 that we use when evaluating later)
                    model_loss(m_idx, v_idx) = -1*mean(R2_val, 2);
                end
                
                % Display model estimation status for user
                if verbose
                    timeNow = datetime('now','TimeZone','local',...
                        'Format','d-MMM-y HH:mm:ss Z');
                    disp(['Done with model order ', num2str(mOrder), ...
                        ' and VNS delay ', num2str(vDelay), ...
                        ' at ', datestr(timeNow)])
                end
            end
        end
    else
        % We can proceed with the usual optimization routine
        % Initialize array to store AICc for each model
        model_loss = zeros(length(modelOrders), length(vnsDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                % Model order to use in this iteration
                mOrder = modelOrders(m_idx);
                
                % VNS delay to use in this iteration
                vDelay =  vnsDelays(v_idx);
                
                % Estimate models using ssest()
                % N4SID estimation followed by prediction error minimization
                modelFit = ssest(fitData, mOrder, ssOptions, ...
                    'InputDelay', vDelay, 'Form', 'modal', 'Ts', sampleTime);
                
                % If the covariance matrix estimate is unreliable or
                % the number of free parameters is less than the number
                % of datapoints available, then I set the AICc value to
                % infinity so that model does not get selected
                if isempty(modelFit.covar) || ...
                        isempty(modelFit.Report.Fit.AICc)
                    model_loss(m_idx, v_idx) = inf;
                else
                    % Store AICc
                    model_loss(m_idx, v_idx) = ...
                        modelFit.Report.Fit.AICc;
                end
                
                % Display model estimation status for user
                if verbose
                    timeNow = datetime('now','TimeZone','local',...
                        'Format','d-MMM-y HH:mm:ss Z');
                    disp(['Done with model order ', num2str(mOrder), ...
                        ' and VNS delay ', num2str(vDelay), ...
                        ' at ', datestr(timeNow)])
                end
            end
        end
        
    end
    
    % Find the optimal model order and delay(s) using my own minimization
    % routine (the min() function in matlab is a bit complicated)
    
    % Initialize
    optLoss = model_loss(1, 1);
    opt_m_idx = 1;
    opt_v_idx = 1;
    
    % Search
    for m_idx = 1:size(model_loss, 1)
        for v_idx = 1:length(vnsDelays)
            % If this AICc is smaller, re-assign
            if model_loss(m_idx, v_idx) < optLoss
                optLoss = model_loss(m_idx, v_idx);
                opt_m_idx = m_idx;
                opt_v_idx = v_idx;
            end
        end
    end
    
    % Store optimal model order
    if max_mOrder
        % Ignore indexing because it does not apply; we only fit one model
        % order: the maximum possible
        mOrder_star = modelOrders(stop_m_idx);
    else
        mOrder_star = modelOrders(opt_m_idx);
    end
    
    % Store delay(s)
    varargout{1} = vnsDelays(opt_v_idx);
    
elseif strcmp(inputOption, 'vt')
    varargout = cell(1, 2);
    vnsDelays = varargin{1};
    traumDelays = varargin{2};
    
    % Check if we're supposed to just use the maximum available
    if max_mOrder
        % Initialize array to store AICc, but only allow one option for
        % model order
        model_loss = zeros(1, length(vnsDelays), length(traumDelays));
        
        % Since we need to use the maximum available model order, we should
        % loop backwards until we fit the model error free
        m_idx = length(modelOrders);    % Initialize at max
        fitSuccess = false;             % Initialize boolean to false
        while ~fitSuccess
            % Start with true for looping through delays to test out
            loopDelays = true; stop_m_idx = m_idx;
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    if loopDelays
                        % Since input delay should have no effect on whether the
                        % model was successfully fit (should have no effect on
                        % number of parameters, which means no consideration of
                        % number of datapoints vs. number of parameters), once the
                        % model was successfully fit once, then we can assume it
                        % will be for all possible delays, so we set fitSuccess to
                        % true, never check again, and start storing away AICc
                        
                        % Model order to use in this iteration
                        mOrder = modelOrders(m_idx);
                        
                        % VNS delay to use in this iteration
                        vDelay =  vnsDelays(v_idx);
                        
                        % Trauma delay to use in this iteration
                        tDelay = traumDelays(t_idx);
                        
                        % Estimate models using ssest()
                        % N4SID estimation followed by prediction error minimization
                        modelFit = ssest(fitData, mOrder, ssOptions, ...
                            'InputDelay', [vDelay; tDelay], 'Form', ...
                            'modal', 'Ts', sampleTime);
                        
                        % If the number of free parameters is less than the number
                        % of datapoints available, then I skip this model order
                        if isempty(modelFit.Report.Fit.AICc)
                            % Set looping boolean to false so we skip the rest of
                            % the delays
                            loopDelays = false;
                        else
                            % If we were successful, then it's time to just stick
                            % with this model order and loop through delays
                            fitSuccess = true;
                        end
                        
                        % If the covariance matrix estimate is unreliable, then I
                        % set the AICc value to infinity so that model does not get
                        % selected (only necessary if loopDelays is still true)
                        if loopDelays && isempty(modelFit.covar)
                            model_loss(1, v_idx, t_idx) = inf;
                        elseif loopDelays
                            % Store AICc
                            model_loss(1, v_idx, t_idx) = ...
                                modelFit.Report.Fit.AICc;
                        end
                        
                        % Display model estimation status for user
                        if verbose
                            timeNow = datetime('now','TimeZone','local',...
                                'Format','d-MMM-y HH:mm:ss Z');
                            disp(['Max model order mode: Done with model order ', ...
                                num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and trauma delay ', num2str(tDelay), ...
                            ' at ', datestr(timeNow)])
                        end
                    end
                end
            end
            % Negative increment
            m_idx = m_idx - 1;
        end
        
    % Check if we need to purposefully leak data for HP optimization
    elseif dataLeak
        % Initialize array to store prediction errors on test data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Trauma delay to use in this iteration
                    tDelay = traumDelays(t_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; tDelay], 'Form', ...
                        'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the loss value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, t_idx) = inf;
                    else
                        % Evaluate model prediction on test data
                        R2_test = eval_modelPred(modelFit, testData, ...
                            'K_steps', 1, 'Rsquared', ...
                            'IC_training', inputArray, featArray);
                        
                        % Store loss on test data (loss will just be negative
                        % of the mean R^2 that we use when evaluating later)
                        model_loss(m_idx, v_idx, t_idx) = -1*mean(R2_test, 2);
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and trauma delay ', num2str(tDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end
        
    % Check if we need to use validation data for HP optimization
    elseif valApproach
        % Initialize array to store prediction errors on val. data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Trauma delay to use in this iteration
                    tDelay = traumDelays(t_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; tDelay], 'Form', ...
                        'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the loss value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, t_idx) = inf;
                    else
                        % Evaluate model prediction on validation data
                        R2_val = eval_modelPred(modelFit, valData, ...
                            'K_steps', 1, 'Rsquared', ...
                            'IC_training', inputArray, featArray);
                        
                        % Store loss on test data (loss will just be negative
                        % of the mean R^2 that we use when evaluating later)
                        model_loss(m_idx, v_idx, t_idx) = -1*mean(R2_val, 2);
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and trauma delay ', num2str(tDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end    
    else
        % We can proceed with the usual optimization routine
        % Initialize array to store AICc for each model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Trauma delay to use in this iteration
                    tDelay = traumDelays(t_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; tDelay], ...
                        'Form', 'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the AICc value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, t_idx) = inf;
                    else
                        % Store AICc
                        model_loss(m_idx, v_idx, t_idx) = ...
                            modelFit.Report.Fit.AICc;
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and trauma delay ', num2str(tDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end
        
    end
    
    % Find the optimal model order and delay(s) using my own minimization
    % routine (the min() function in matlab is a bit complicated)
    
    % Initialize
    optLoss = model_loss(1, 1, 1);
    opt_m_idx = 1;
    opt_v_idx = 1;
    opt_t_idx = 1;
    
    % Search
    for m_idx = 1:size(model_loss, 1)
        for v_idx = 1:length(vnsDelays)
            for t_idx = 1:length(traumDelays)
                % If this AICc is smaller, re-assign
                if model_loss(m_idx, v_idx, t_idx) < optLoss
                    optLoss = model_loss(m_idx, v_idx, t_idx);
                    opt_m_idx = m_idx;
                    opt_v_idx = v_idx;
                    opt_t_idx = t_idx;
                end
            end
        end
    end
    
    % Store optimal model order
    if max_mOrder
        % Ignore indexing because it does not apply; we only fit one model
        % order: the maximum possible
        mOrder_star = modelOrders(stop_m_idx);
    else
        mOrder_star = modelOrders(opt_m_idx);
    end
    
    % Store delays
    varargout{1} = vnsDelays(opt_v_idx);
    varargout{2} = traumDelays(opt_t_idx);
    
elseif strcmp(inputOption, 'vn')
    varargout = cell(1, 2);
    vnsDelays = varargin{1};
    neutrDelays = varargin{2};
    
    % Check if we're supposed to just use the maximum available
    if max_mOrder
        % Initialize array to store AICc, but only allow one option for
        % model order
        model_loss = zeros(1, length(vnsDelays), length(neutrDelays));
        
        % Since we need to use the maximum available model order, we should
        % loop backwards until we fit the model error free
        m_idx = length(modelOrders);    % Initialize at max
        fitSuccess = false;             % Initialize boolean to false
        while ~fitSuccess
            % Start with true for looping through delays to test out
            loopDelays = true; stop_m_idx = m_idx;
            for v_idx = 1:length(vnsDelays)
                for n_idx = 1:length(neutrDelays)
                    if loopDelays
                        % Since input delay should have no effect on whether the
                        % model was successfully fit (should have no effect on
                        % number of parameters, which means no consideration of
                        % number of datapoints vs. number of parameters), once the
                        % model was successfully fit once, then we can assume it
                        % will be for all possible delays, so we set fitSuccess to
                        % true, never check again, and start storing away AICc
                        
                        % Model order to use in this iteration
                        mOrder = modelOrders(m_idx);
                        
                        % VNS delay to use in this iteration
                        vDelay =  vnsDelays(v_idx);
                        
                        % Trauma delay to use in this iteration
                        nDelay = neutrDelays(n_idx);
                        
                        % Estimate models using ssest()
                        % N4SID estimation followed by prediction error minimization
                        modelFit = ssest(fitData, mOrder, ssOptions, ...
                            'InputDelay', [vDelay; nDelay], 'Form', ...
                            'modal', 'Ts', sampleTime);
                        
                        % If the number of free parameters is less than the number
                        % of datapoints available, then I skip this model order
                        if isempty(modelFit.Report.Fit.AICc)
                            % Set looping boolean to false so we skip the rest of
                            % the delays
                            loopDelays = false;
                        else
                            % If we were successful, then it's time to just stick
                            % with this model order and loop through delays
                            fitSuccess = true;
                        end
                        
                        % If the covariance matrix estimate is unreliable, then I
                        % set the AICc value to infinity so that model does not get
                        % selected (only necessary if loopDelays is still true)
                        if loopDelays && isempty(modelFit.covar)
                            model_loss(1, v_idx, n_idx) = inf;
                        elseif loopDelays
                            % Store AICc
                            model_loss(1, v_idx, n_idx) = ...
                                modelFit.Report.Fit.AICc;
                        end
                        
                        % Display model estimation status for user
                        if verbose
                            timeNow = datetime('now','TimeZone','local',...
                                'Format','d-MMM-y HH:mm:ss Z');
                            disp(['Max model order mode: Done with model order ', ...
                                num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and neutral delay ', num2str(nDelay), ...
                            ' at ', datestr(timeNow)])
                        end
                    end
                end
            end
            % Negative increment
            m_idx = m_idx - 1;
        end
    
        
    % Check if we need to purposefully leak data for HP optimization
    elseif dataLeak
        % Initialize array to store prediction errors on test data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for n_idx = 1:length(neutrDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Trauma delay to use in this iteration
                    nDelay = neutrDelays(n_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; nDelay], 'Form', ...
                        'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the loss value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, n_idx) = inf;
                    else
                        % Evaluate model prediction on test data
                        R2_test = eval_modelPred(modelFit, testData, ...
                            'K_steps', 1, 'Rsquared', ...
                            'IC_training', inputArray, featArray);
                        
                        % Store loss on test data (loss will just be negative
                        % of the mean R^2 that we use when evaluating later)
                        model_loss(m_idx, v_idx, n_idx) = -1*mean(R2_test, 2);
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and neutral delay ', num2str(tDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end
    
    % Check if we need to use validation data for HP optimization
    elseif valApproach
        % Initialize array to store prediction errors on val. data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for n_idx = 1:length(neutrDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Trauma delay to use in this iteration
                    nDelay = neutrDelays(n_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; nDelay], 'Form', ...
                        'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the loss value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, n_idx) = inf;
                    else
                        % Evaluate model prediction on validation data
                        R2_val = eval_modelPred(modelFit, valData, ...
                            'K_steps', 1, 'Rsquared', ...
                            'IC_training', inputArray, featArray);
                        
                        % Store loss on test data (loss will just be negative
                        % of the mean R^2 that we use when evaluating later)
                        model_loss(m_idx, v_idx, n_idx) = -1*mean(R2_val, 2);
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and neutral delay ', num2str(nDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end    
    else
        % We can proceed with the usual optimization routine
        % Initialize array to store AICc for each model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for n_idx = 1:length(neutrDelays)
                    % Model order to use in this iteration
                    mOrder = modelOrders(m_idx);
                    
                    % VNS delay to use in this iteration
                    vDelay =  vnsDelays(v_idx);
                    
                    % Neutral delay to use in this iteration
                    nDelay = neutrDelays(n_idx);
                    
                    % Estimate models using ssest()
                    % N4SID estimation followed by prediction error minimization
                    modelFit = ssest(fitData, mOrder, ssOptions, ...
                        'InputDelay', [vDelay; nDelay], ...
                        'Form', 'modal', 'Ts', sampleTime);
                    
                    % If the covariance matrix estimate is unreliable or
                    % the number of free parameters is less than the number
                    % of datapoints available, then I set the AICc value to
                    % infinity so that model does not get selected
                    if isempty(modelFit.covar) || ...
                            isempty(modelFit.Report.Fit.AICc)
                        model_loss(m_idx, v_idx, n_idx) = inf;
                    else
                        % Store AICc
                        model_loss(m_idx, v_idx, n_idx) = ...
                            modelFit.Report.Fit.AICc;
                    end
                    
                    % Display model estimation status for user
                    if verbose
                        timeNow = datetime('now','TimeZone','local',...
                            'Format','d-MMM-y HH:mm:ss Z');
                        disp(['Done with model order ', num2str(mOrder), ...
                            ' and VNS delay ', num2str(vDelay), ...
                            ' and neutral delay ', num2str(nDelay), ...
                            ' at ', datestr(timeNow)])
                    end
                end
            end
        end
        
    end
    
    % Find the optimal model order and delay(s) using my own minimization
    % routine (the min() function in matlab is a bit complicated)
    
    % Initialize
    optLoss = model_loss(1, 1, 1);
    opt_m_idx = 1;
    opt_v_idx = 1;
    opt_n_idx = 1;
    
    % Search
    for m_idx = 1:size(model_loss, 1)
        for v_idx = 1:length(vnsDelays)
            for n_idx = 1:length(neutrDelays)
                % If this AICc is smaller, re-assign
                if model_loss(m_idx, v_idx, n_idx) < optLoss
                    optLoss = model_loss(m_idx, v_idx, n_idx);
                    opt_m_idx = m_idx;
                    opt_v_idx = v_idx;
                    opt_n_idx = n_idx;
                end
            end
        end
    end
    
    % Store optimal model order
    if max_mOrder
        % Ignore indexing because it does not apply; we only fit one model
        % order: the maximum possible
        mOrder_star = modelOrders(stop_m_idx);
    else
        mOrder_star = modelOrders(opt_m_idx);
    end
    
    % Store delays
    varargout{1} = vnsDelays(opt_v_idx);
    varargout{2} = neutrDelays(opt_n_idx);
    
else
    varargout = cell(1, 3);
    vnsDelays = varargin{1};
    traumDelays = varargin{2};
    neutrDelays = varargin{3};
    
    % Check if we're supposed to just use the maximum available
    if max_mOrder
        % Initialize array to store AICc, but only allow one option for
        % model order
        model_loss = zeros(1, length(vnsDelays), length(traumDelays), ...
            length(neutrDelays));
        
        % Since we need to use the maximum available model order, we should
        % loop backwards until we fit the model error free
        m_idx = length(modelOrders);    % Initialize at max
        fitSuccess = false;             % Initialize boolean to false
        while ~fitSuccess
            % Start with true for looping through delays to test out
            loopDelays = true; stop_m_idx = m_idx;
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    for n_idx = 1:length(neutrDelays)
                        if loopDelays
                            % Since input delay should have no effect on whether the
                            % model was successfully fit (should have no effect on
                            % number of parameters, which means no consideration of
                            % number of datapoints vs. number of parameters), once the
                            % model was successfully fit once, then we can assume it
                            % will be for all possible delays, so we set fitSuccess to
                            % true, never check again, and start storing away AICc
                            
                            % Model order to use in this iteration
                            mOrder = modelOrders(m_idx);
                            
                            % VNS delay to use in this iteration
                            vDelay =  vnsDelays(v_idx);
                            
                            % Trauma delay to use in this iteration
                            tDelay = traumDelays(t_idx);
                            
                            % Neutral delay to use in this iteration
                            nDelay = neutrDelays(n_idx);
                            
                            % Estimate models using ssest()
                            % N4SID estimation followed by prediction error minimization
                            modelFit = ssest(fitData, mOrder, ssOptions, ...
                                'InputDelay', [vDelay; tDelay; nDelay], 'Form', ...
                                'modal', 'Ts', sampleTime);
                            
                            % If the number of free parameters is less than the number
                            % of datapoints available, then I skip this model order
                            if isempty(modelFit.Report.Fit.AICc)
                                % Set looping boolean to false so we skip the rest of
                                % the delays
                                loopDelays = false;
                            else
                                % If we were successful, then it's time to just stick
                                % with this model order and loop through delays
                                fitSuccess = true;
                            end
                            
                            % If the covariance matrix estimate is unreliable, then I
                            % set the AICc value to infinity so that model does not get
                            % selected (only necessary if loopDelays is still true)
                            if loopDelays && isempty(modelFit.covar)
                                model_loss(1, v_idx, t_idx, n_idx) = inf;
                            elseif loopDelays
                                % Store AICc
                                model_loss(1, v_idx, t_idx, n_idx) = ...
                                    modelFit.Report.Fit.AICc;
                            end
                            
                            % Display model estimation status for user
                            if verbose
                                timeNow = datetime('now','TimeZone','local',...
                                    'Format','d-MMM-y HH:mm:ss Z');
                                disp(['Max model order mode: Done with model order ', ...
                                    num2str(mOrder), ...
                                    ' and VNS delay ', num2str(vDelay), ...
                                    ' and trauma delay ', num2str(tDelay), ...
                                    ' and neutral delay ', num2str(nDelay), ...
                                    ' at ', datestr(timeNow)])
                            end
                        end
                    end
                    
                end
            end
            % Negative increment
            m_idx = m_idx - 1;
        end
        
    % Check if we need to purposefully leak data for HP optimization
    elseif dataLeak
        % Initialize array to store prediction errors on test data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    for n_idx = 1:length(neutrDelays)
                        % Model order to use in this iteration
                        mOrder = modelOrders(m_idx);
                        
                        % VNS delay to use in this iteration
                        vDelay =  vnsDelays(v_idx);
                        
                        % Trauma delay to use in this iteration
                        tDelay = traumDelays(t_idx);
                        
                        % Trauma delay to use in this iteration
                        nDelay = traumDelays(n_idx);
                        
                        % Estimate models using ssest()
                        % N4SID estimation followed by prediction error minimization
                        modelFit = ssest(fitData, mOrder, ssOptions, ...
                            'InputDelay', [vDelay; tDelay; nDelay], 'Form', ...
                            'modal', 'Ts', sampleTime);
                        
                        % If the covariance matrix estimate is unreliable or
                        % the number of free parameters is less than the number
                        % of datapoints available, then I set the loss value to
                        % infinity so that model does not get selected
                        if isempty(modelFit.covar) || ...
                                isempty(modelFit.Report.Fit.AICc)
                            model_loss(m_idx, v_idx, t_idx, n_idx) = inf;
                        else
                            % Evaluate model prediction on test data
                            R2_test = eval_modelPred(modelFit, testData, ...
                                'K_steps', 1, 'Rsquared', ...
                                'IC_training', inputArray, featArray);
                            
                            % Store loss on test data (loss will just be negative
                            % of the mean R^2 that we use when evaluating later)
                            model_loss(m_idx, v_idx, t_idx, n_idx) = ...
                                -1*mean(R2_test, 2);
                        end
                        
                        % Display model estimation status for user
                        if verbose
                            timeNow = datetime('now','TimeZone','local',...
                                'Format','d-MMM-y HH:mm:ss Z');
                            disp(['Done with model order ', num2str(mOrder), ...
                                ' and VNS delay ', num2str(vDelay), ...
                                ' and trauma delay ', num2str(tDelay), ...
                                ' and neutral delay ', num2str(tDelay), ...
                                ' at ', datestr(timeNow)])
                        end
                    end
                end
            end
        end
    
    % Check if we need to use validation data for HP optimization
    elseif valApproach
        % Initialize array to store prediction errors on val. data for each
        % model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    for n_idx = 1:length(neutrDelays)
                        % Model order to use in this iteration
                        mOrder = modelOrders(m_idx);
                        
                        % VNS delay to use in this iteration
                        vDelay =  vnsDelays(v_idx);
                        
                        % Trauma delay to use in this iteration
                        tDelay = traumDelays(t_idx);
                        
                        % Trauma delay to use in this iteration
                        nDelay = traumDelays(n_idx);
                        
                        % Estimate models using ssest()
                        % N4SID estimation followed by prediction error minimization
                        modelFit = ssest(fitData, mOrder, ssOptions, ...
                            'InputDelay', [vDelay; tDelay; nDelay], 'Form', ...
                            'modal', 'Ts', sampleTime);
                        
                        % If the covariance matrix estimate is unreliable or
                        % the number of free parameters is less than the number
                        % of datapoints available, then I set the loss value to
                        % infinity so that model does not get selected
                        if isempty(modelFit.covar) || ...
                                isempty(modelFit.Report.Fit.AICc)
                            model_loss(m_idx, v_idx, t_idx, n_idx) = inf;
                        else
                            % Evaluate model prediction on validation data
                            R2_val = eval_modelPred(modelFit, valData, ...
                                'K_steps', 1, 'Rsquared', ...
                                'IC_training', inputArray, featArray);
                            
                            % Store loss on test data (loss will just be negative
                            % of the mean R^2 that we use when evaluating later)
                            model_loss(m_idx, v_idx, t_idx, n_idx) = ...
                                -1*mean(R2_val, 2);
                        end
                        
                        % Display model estimation status for user
                        if verbose
                            timeNow = datetime('now','TimeZone','local',...
                                'Format','d-MMM-y HH:mm:ss Z');
                            disp(['Done with model order ', num2str(mOrder), ...
                                ' and VNS delay ', num2str(vDelay), ...
                                ' and trauma delay ', num2str(tDelay), ...
                                ' and neutral delay ', num2str(tDelay), ...
                                ' at ', datestr(timeNow)])
                        end
                    end
                end
            end
        end
        
    else
        % We can proceed with the usual optimization routine
        % Initialize array to store AICc for each model
        model_loss = zeros(length(modelOrders), ...
            length(vnsDelays), length(traumDelays), length(neutrDelays));
        
        % Loop through model orders and delays to optimize over
        for m_idx = 1:length(modelOrders)
            for v_idx = 1:length(vnsDelays)
                for t_idx = 1:length(traumDelays)
                    for n_idx = 1:length(neutrDelays)
                        % Model order to use in this iteration
                        mOrder = modelOrders(m_idx);
                        
                        % VNS delay to use in this iteration
                        vDelay =  vnsDelays(v_idx);
                        
                        % Trauma delay to use in this iteration
                        tDelay = traumDelays(t_idx);
                        
                        % Neutral dleay to use in this iteration
                        nDelay = neutrDelays(n_idx);
                        
                        % Estimate models using ssest()
                        % N4SID estimation followed by prediction error min.
                        modelFit = ssest(fitData, mOrder, ssOptions, ...
                            'InputDelay', [vDelay; tDelay; nDelay], ...
                            'Form', 'modal', 'Ts', sampleTime);
                        
                        
                        % If the covariance matrix estimate is unreliable or
                        % the number of free parameters is less than the number
                        % of datapoints available, then I set the AICc value to
                        % infinity so that model does not get selected
                        if isempty(modelFit.covar) || ...
                                isempty(modelFit.Report.Fit.AICc)
                            model_loss(m_idx, v_idx, t_idx, n_idx) = inf;
                        else
                            % Store AICc
                            model_loss(m_idx, v_idx, t_idx, n_idx) = ...
                                modelFit.Report.Fit.AICc;
                        end
                        
                        % Display model estimation status for user
                        if verbose
                            timeNow = datetime('now','TimeZone','local',...
                                'Format','d-MMM-y HH:mm:ss Z');
                            disp(['Done with model order ', num2str(mOrder), ...
                                ' and VNS delay ', num2str(vDelay), ...
                                ' and trauma delay ', num2str(tDelay), ...
                                ' and neutral delay ', num2str(nDelay), ...
                                ' at ', datestr(timeNow)])
                        end
                    end
                end
            end
        end
        
    end
    
    % Find the optimal model order and delay(s) using my own minimization
    % routine (the min() function in matlab is a bit complicated)
    
    % Initialize
    optLoss = model_loss(1, 1, 1, 1);
    opt_m_idx = 1;
    opt_v_idx = 1;
    opt_t_idx = 1;
    opt_n_idx = 1;
    
    % Search
    for m_idx = 1:size(model_loss, 1)
        for v_idx = 1:length(vnsDelays)
            for t_idx = 1:length(traumDelays)
                for n_idx = 1:length(neutrDelays)
                    % If this AICc is smaller, re-assign
                    if model_loss(m_idx, v_idx, t_idx, n_idx) < optLoss
                        optLoss = model_loss(m_idx, v_idx, t_idx, n_idx);
                        opt_m_idx = m_idx;
                        opt_v_idx = v_idx;
                        opt_t_idx = t_idx;
                        opt_n_idx = n_idx;
                    end
                end
            end
        end
    end
    
    % Store optimal model order
    if max_mOrder
        % Ignore indexing because it does not apply; we only fit one model
        % order: the maximum possible
        mOrder_star = modelOrders(stop_m_idx);
    else
        mOrder_star = modelOrders(opt_m_idx);
    end
    
    % Store delay(s)
    varargout{1} = vnsDelays(opt_v_idx);
    varargout{2} = traumDelays(opt_t_idx);
    varargout{3} = neutrDelays(opt_n_idx);
end

end

