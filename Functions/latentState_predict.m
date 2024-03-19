function xPredict = ...
    latentState_predict(ssModel, yData, uData, Ts, varargin)
% Because the built-in predict() function provided by MATLAB for identified
% state-space models only returns predicted outputs, y, and not the
% predicted latent state, x, this function uses the workaround recommended
% by MATLAB to produce predicted response using simulation. The key here
% being that the simulation function, sim(), returns the latent state, x

% ------------ Inputs --------------- %
% ssModel: idss object that defines state-space model
% yData: output measurement array of TxM dimension, where T is number of 
%        timesteps and M is number out output channels
% uData: input array of TxU dimension, where T is number of timesteps and 
%        U is number of input channels; if autonomous system, enter empty 
%        array (i.e., [])
% Ts: sample time

% (optional)
% 'K_steps': flag followed by number of steps ahead to predict (default = 1)
% 'predictOpt': flag followed by predictOptions() object


% ----------- Outputs ---------------- %
% xPredict: TxN array of latent state predictions, where T is number of
%           timesteps and N is dimensionality of latent state


% Parse varargin
for arg = 1:length(varargin)
    if strcmp(varargin{arg}, 'K_steps')
        K_steps = varargin{arg + 1};
    elseif strcmp(varargin{arg}, 'predictOpt')
        predictOpt = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('K_steps', 'var')
    K_steps = 1;
end
if ~exist('predictOpt', 'var')
    predictOpt = predictOptions;
end

% First thing to do is to create an iddata object from yData, uData, and Ts
dataObj = iddata(yData, uData, Ts);

% Now we can apply the predict() function to obtain our initial conditions
% and predictor model to use during simulation (I also return y_predict to
% do a sanity check later)
[y_predict, initialCond, predModel] = predict(ssModel, dataObj, ...
    K_steps, predictOpt);

% Specify the returned initial conditions as the initial conditions to use
% during simulation
simOpt = simOptions('InitialCondition', initialCond);

% Simulate the predictor model forward as done in documentation
% (https://www.mathworks.com/help/ident/ref/predict.html)
[y_sim, ~, xPredict, ~] = sim(predModel, ...
    [dataObj.OutputData, dataObj.InputData], simOpt);



% ---------- Sanity checks ----------- %
% Compute difference between the two arrays that are supposed to be equal
y_residual = y_predict.OutputData - y_sim;

% Sum up difference (quick and dirty way for comparison)
sumY = sum(y_residual, 'all');

% Make sure sum of differences is 0, i.e., negligible compared to smallest
% magnitudes relevant to data
if sumY > min(abs(y_sim), [], 'all')/1000
    error('Something went wrong with y')
end


% Transpose initial conditions if needed for comparison
if size(initialCond, 1) > 1
    initialCond = initialCond';
end

% Compute difference
x0_residual = xPredict(1, :) - initialCond;

% Sum up difference (quick and dirty way for comparison)
sumX0 = sum(x0_residual);

% Make sure sum of differences is 0, i.e., negligible compared to smallest
% magnitudes relevant to data
if sumX0 > min(abs(xPredict), [], 'all')/1000
    error('Something went wrong with x_0')
end


end

