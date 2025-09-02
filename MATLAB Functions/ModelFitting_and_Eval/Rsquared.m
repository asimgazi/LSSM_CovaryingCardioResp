function R2 = Rsquared(actual, predicted)
% This function takes a vector or array of actual values and a vector or 
% array of predicted values and computes the coefficient of determination
% per pair of columns (single R^2 will be returned if both are vectors

%------------ Inputs ------------- %
% actual: If vector, Nx1 or 1xN vector of actual values
%         If array, NxK array of actual values, where K is the number of
%         outputs
% predicted: If vector, Nx1 or 1xN vector of predicted values
%            If array, NxK array of predicted values, where K is the number
%            of outputs
% Note that we must have N > 1

% ------- Outputs ---------- %
% R2: If actual and predicted were vectors, this will be a scalar R^2 
%     value, computed using the formula 1 - MSE/Variance, where variance is
%     computed by normalizing by N, rather than N-1
%     If arrays were input (i.e., multiple outputs were provided), a 1xK
%     vector of R^2 will be returned

% Check both dimensions of actual and predicted; if they're vectors, make
% them both column vectors so below code can be the same
dim1_a = size(actual, 1); dim1_p = size(predicted, 1);
if dim1_a == 1
    actual = actual';
end
if dim1_p == 1
    predicted = predicted';
end

% Initialize R2
R2 = zeros(1, size(actual, 2));

% Loop through columns, computing R2 one at a time
for r = 1:length(R2)
    % Let's take care of the easy part of calculating the data's variance,
    % normalized by N instead of N-1 so the Ns cancel out from the mean square
    % error
    
    % Inputting a 1 for weight, w, specifies normalization by N
    denom = var(actual(:, r), 1);
    
    % Now let's compute the mean square error
    numer = mean((predicted(:, r) - actual(:, r)).^2);
    
    % Compute the final R^2 value
    R2(r) = 1 - numer/denom;
end

end

