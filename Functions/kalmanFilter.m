function [x_hat, y_hat, varargout] = kalmanFilter(A, C, Q, R, S, y, varargin)
% This function takes system matrices A, C and covariance matrices Q, R, S,
% and performs Kalman filtering on the observations y, to produce state
% estimates at each corresponding time point, x_hat
% This function can be used for as many steps as desired, as it simply
% iterates through the rows present in y

% Inputs
% A: state transition matrix - N x N
% C: state to output mapping matrix - M x N
% Q: process noise covariance - N x N
% R: measurement noise covariance - M x M
% S: cross-covariance between two noises - N x M
% y: output sequence - T x M, where T is number of time instances
% (optional)
% 'forLoop' - flag to specify whether for loop approach is being used
% x_pred_in: predicted state from previous step
% P_pred_in: predicted error covariance from previous step
% 'input' - flag to specify whether an input to the system exists
% B: input to state mapping matrix - N x K
% u: input sequence - T x K, where T is number of time instances

% Outputs
% x_hat: estimated state sequence - T x N
% y_hat: estimated output sequence - T x M
% (optional)
% x_pred_out: predicted state for next step (varargout{1})
% y_pred_out: predicted output for next step (varargout{2})
% P_pred_out: predicted error covariance for next step (varargout{3})

% Initialize defaults
forLoop = false;
nonAut = false;

% Parse varargin
if ~isempty(varargin)
    for arg = 1:length(varargin)
        if strcmp(varargin{arg}, 'forLoop')
            forLoop = true;
            x_pred_in = varargin{arg + 1};
            P_pred_in = varargin{arg + 2};
        elseif strcmp(varargin{arg}, 'input')
            nonAut = true;  % Non-autonomous case
            B = varargin{arg + 1};
            u = varargin{arg + 2};
        end
    end
end

% System Equations (for reference)
% x_i+1 = A*x_i + B*u_i + w_i
% y_i = C*x_i + v_i
% E[w_i*w_i'] = Q
% E[v_i*v_i'] = R
% E[w_i*v_i'] = S

% Herein, we assume S = 0
if S(1) ~= 0
    disp('Function was not made for correlated noises')
end


% Initialize prediction and measurement update state vectors
x_hat = zeros(size(y, 1), size(A, 1));
x_pred = zeros(size(y, 1), size(A, 1));

% Initialize error covariance matrix
P_hat = zeros(size(y, 1), size(A, 1), size(A, 1));
P_pred = zeros(size(y, 1), size(A, 1), size(A, 1));


if forLoop
    % Initialize first prediction to the inputs to the function
    x_pred(1, :) = x_pred_in;
    P_pred(1, :, :) = P_pred_in;
else
    % Initialize first prediction to first observation, just to initialize
    % (essentially just implement naive prediction for first time step)
    x_pred(1, :) = y(1, :);
    P_pred(1, :, :) = zeros(size(A, 1), size(A, 1));
end


% For every timestep provided
for t = 1:size(y, 1)
    % Measurement update
    % (note that transposes had to be used in certain places to change the
    % row instances into column vectors)
    Kalm_gain = squeeze(P_pred(t, :, :))*C'*inv(C*squeeze(P_pred(t, :, :))*C' + R);
    x_hat(t, :) = (x_pred(t, :)' + Kalm_gain*(y(t, :)' - C*x_pred(t, :)'))';
    P_hat(t, :, :) = (eye(size(A, 1)) - Kalm_gain*C)*squeeze(P_pred(t, :, :));
    
    % Prediction update
    x_pred(t+1, :) = (A*x_hat(t, :)')';
    P_pred(t+1, :, :) = A*squeeze(P_hat(t, :, :))*A' + Q;
end

% The final prediction step actually appends to x_pred and P_pred beyond
% their original initializations. Hence, we can use that to our advantage
% if we would like to store that final prediction
if forLoop
    x_pred_out = x_pred(end, :);
    P_pred_out = squeeze(P_pred(end, :, :));
    
    varargout{1} = x_pred_out;
    varargout{2} = P_pred_out;  
end

end

