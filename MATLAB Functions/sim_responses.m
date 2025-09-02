function [outputs, states] = sim_responses(ssModel, input_sim, inputsOpt)
% This function takes a state-space model, an input signal for simulation
% and the option for protocol condition inputs and simulates forward
% outputs and states for all conditions treated as inputs

% ------------- Inputs -------------- %
% ssModel: idss object that stores state-space model to simulate forward
% input_sim: Either Tx1 array that stores input to use for each condition
%            separately, OR Tx? array where an input signal to use for
%            simulation is provided for each protocol condition used as
%            input (? because number depends on inputsOpt); NOTE THAT THERE
%            IS AN ADDITIONAL IMPLICATION OF 1 VS. ? DIMENSION: *R*
%            (defined below)
%            If you want to just simulate a single response (R = 1), then 
%            you need to specify the input for each condition. If you want 
%            to simulate multiple reponses, then you need to input a Tx1
% inputsOpt: string that specifies which protocol conditions
%            inputs were created for - 
%            'v': just VNS
%            'vt': VNS and trauma recall
%            'vn': VNS and neutral condition
%            'vtn': VNS, trauma recall, and neutral

% ------------ Outputs -------------- %
% outputs: TxMxR array, where M is the number of features, T is the number
%          of timepoints, and R is the number of simulated repsonses. This
%          array stores the outputs of the state-space model (often denoted
%          y)
% states: TxNxR array. where N is the number of states, T is the number of
%         timepoints, and R is the number of simulated responses. This
%         array stores the states of the state-space model (often denoted
%         x)


% Store number of inputs based on inputsOpt
if strcmp(inputsOpt, 'vtn')
    numInputs = 3;
elseif strcmp(inputsOpt, 'vt') || strcmp(inputsOpt, 'vn')
    numInputs = 2;
else
    numInputs = 1;
end

% Is our input of the required dimension?
if size(input_sim, 2) < numInputs
    % Create multiple inputs of appropriate dimension (MATLAB System ID
    % toolbox expects each input channel to be its own *column*, not row
    inputSet = zeros(size(input_sim, 1), numInputs, numInputs);
    
    % Loop through all responses required
    for r = 1:size(inputSet, 3)
        % Store simulation input for appropriate condition's column
        inputSet(:, r, r) = input_sim;
    end
else
    % Store the single input (again, MATLAB System ID toolbox expects each
    % input channel to be its own *column*, not row)
    inputSet = input_sim;
end

% Initialize and then later concatenate
outputs = []; states = [];

% Loop through all responses required
for r = 1:size(inputSet, 3)
    % The ignored outputs are y_sd and x_sd, i.e., the point by point 
    % standard deviation estimates of each element of the output and state,
    % respectively
    [y, ~, x, ~] = sim(ssModel, inputSet(:, :, r));
    
    % Concatenate along 3rd dimension (depth stores number of responses, R)
    outputs = cat(3, outputs, y); states = cat(3, states, x);
end

end

