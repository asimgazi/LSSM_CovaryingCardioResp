function metric_vec = eval_modelPred(trainedModel, testData, varargin)
% This function takes a model and evaluates its prediction ability on a
% given set of data based on the evaluation options provided

% --------- Inputs ----------- %
% trainedModel: model object that stores the trained model to be evaluated
% testData: iddata object that stores the data to be tested on
% (Optional)
% 'K_steps': flag that specifies that we will be inputting the number of
%            steps to predict ahead for evaluating the model. Next input
%            needs to be the number of steps (e.g., 3 for 3-step ahead
%            prediction). Input inf for infinite step ahead prediction
%            (i.e., stimulation)
% 'IC_training': flag that specifies that initial conditions should be set
%                based on the end of training data. The next two inputs
%                need to be array of inputs (NxK) and array of features in
%                that order (NxM), where K is number of inputs and M is
%                numebr of features
% 'Rsquared': flag that specifies to evaluate the model using coefficient
%             of determination R^2

% -------- Outputs ---------- %
% metric_vec: 1xM vector of metric values, one for each feature

% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'K_steps')
        % Store steps to predict ahead (inf will be passed for simulation)
        K_steps = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'IC_training')
        % Set the initial condition choice accordingly
        IC_choice = 1;
        
        % Store input data and features
        inputTrain = varargin{arg + 1};
        outputTrain = varargin{arg + 2};
    elseif strcmp(varargin{arg}, 'Rsquared')
        % Set the metric choice accordingly
        metric_choice = 0;
    end
end

% Set defaults
if ~exist('K_steps', 'var'); K_steps = 1; end
if ~exist('IC_choice', 'var'); IC_choice = 1; end
if ~exist('metric_choice', 'var'); metric_choice = 0; end

if IC_choice == 1
    % We need to set the initial conditions for prediction based on the end
    % of data provided
    training_struct = struct('Input', inputTrain, 'Output', outputTrain);
    predictOpt = predictOptions('InitialCondition', training_struct);
else
    % Set the initial conditions based on default approach
    predictOpt = predictOptions;
end

% Store away predicted feature values on test set
[y_hat, ~, ~] = predict(trainedModel, testData, K_steps, predictOpt);

if metric_choice == 1
    % Undefined at the moment, but left here in case
else
    % Compute R^2 for each feature
    metric_vec = Rsquared(testData.OutputData, y_hat.OutputData);
end


end

