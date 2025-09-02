function output_struct = segm_condResponses(yArr, tVec, protTable, ...
    t_befStart, t_aftSeg, varargin)
% This function takes an array of feature data, a time vector, and protocol
% timings for a particular subject and segments the data according to
% protocol timings and the time desired before and after each condition

% ------------- Inputs -------------- %
% yArr: TxM array of feature data, where T is number of time points and M
%       is number of features
% tVec: Tx1 vector of time points
% protTable: table of protocol timings
% t_befStart: time to keep before protocol start time when segmenting
% t_aftSeg: time to keep after the segment when segmenting
% (Optional)
% t_vns_seg: time to treat as VNS segment
% t_traum_seg: time to treat as trauma segment
% t_neutr_seg: time to treat as neutral segment
% (total time of a segment will be equal to 
% t_befStart + t_[cond]_seg + t_aftSeg)
% 'sepBaselines' - flag that tells the function to use a specific amount of
%                  time before each segmented response to "renormalize",
%                  i.e., set a separate baseline for each segmented
%                  response to focus on the transient response change
%                  specific to the condition; this flag needs to be
%                  followed by the amount of time to use before each
%                  response to average for the separate baselines

% ------------ Outputs -------------- %
% output_struct: struct with three fields, one for each protocol condition
%               (VNS: v, trauma: t, neutral: n). Each struct will store a
%               NxMxK array, where N is number of time points in a segment,
%               M is number of features, and K is number of condition
%               repetitions (i.e., neutral 2 data will be stored such that
%               dim 3 = 2)

% Since this function is written specifically for DARPA PTSD data, we
% include the protocol ordering as a reminder of what we should expect (and
% if this data are missing, we include NaN)
% N1|N2|T1V1|T2V2|V3|V4|N3|N4|T3V5|T4V6|Lunch
numV = 6;
numT = 4;
numN = 4;

% Parse varargin
for arg = 1:length(varargin)
    if arg == 1
        t_vns_seg = varargin{arg};
    elseif arg == 2
        t_traum_seg = varargin{arg};
    elseif arg == 3
        t_neutr_seg = varargin{arg};
    elseif strcmp(varargin{arg}, 'sepBaselines')
        sepBaselines = true;
        t_sepBaseline = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('t_vns_seg', 'var'); t_vns_seg = 120; end
if ~exist('t_traum_seg', 'var'); t_traum_seg = 60; end
if ~exist('t_neutr_seg', 'var'); t_neutr_seg = 60; end
if ~exist('sepBaselines', 'var'); sepBaselines = false; end
if ~exist('t_sepBaseline', 'var'); t_sepBaseline = 30; end

% Figure out resampling frequency to know how to initialize arrays
Fs = 1/(tVec(2) - tVec(1));

% Initialize the output struct exactly as I want it
% (recall that I will fill missing conditions with NaN, so I know exactly
% how many repetitions I need to initialize)
output_struct = struct('v', ...
    zeros(round((t_befStart + t_vns_seg + t_aftSeg)*Fs), ...
    size(yArr, 2), numV), ...
    't', zeros(round((t_befStart + t_traum_seg + t_aftSeg)*Fs), ...
    size(yArr, 2), numT), ...
    'n', zeros(round((t_befStart + t_neutr_seg + t_aftSeg)*Fs), ...
    size(yArr, 2), numN));

% ---------- Take care of VNS/Sham ------------- %
% Store VNS start and end times for segmentation
vns_start = table2array(protTable(:, 'vns_start')) - t_befStart;
vns_end = vns_start + t_befStart + t_vns_seg + t_aftSeg;
for v_idx = 1:numV
    % If this participant had missing timings, then that means they ended
    % the protocol early; just store NaNs
    if v_idx > length(vns_start)
        output_struct.v(:, :, v_idx) = ...
            nan*ones(size(output_struct.v(:, :, v_idx)));
    else
        % If a NaN is stored, ignore this timing and store NaNs
        if isnan(vns_start(v_idx))
            output_struct.v(:, :, v_idx) = ...
                nan*ones(size(output_struct.v(:, :, v_idx)));
        else
            % Isolate the segment of data corresponding to this timing
            dataSeg = yArr(tVec >= vns_start(v_idx) & ...
                tVec <= vns_end(v_idx), :);
            
            % If we need a separate baseline for each response
            if sepBaselines
                % Compute baseline for this segment of data
                baseline = mean(yArr(...
                    tVec >= vns_start(v_idx) + t_befStart - t_sepBaseline ...
                    & tVec <= vns_start(v_idx) + t_befStart, :), 1);
                
                % Subtract away the baseline
                for feat_idx = 1:size(yArr, 2)
                    dataSeg(:, feat_idx) = dataSeg(:, feat_idx) - ...
                        baseline(feat_idx);
                end
            end
            
            % Make sure enough elements exist in this segment
            if size(dataSeg, 1) == size(output_struct.v, 1)
                % Store away
                output_struct.v(:, :, v_idx) = dataSeg;
            elseif size(dataSeg, 1) == size(output_struct.v, 1) - 1
                % Here, we handle edge case of having one less datapoint
                % We just repeat the final datapoint
                output_struct.v(:, :, v_idx) = [dataSeg; dataSeg(end, :)];
            elseif size(dataSeg, 1) == size(output_struct.v, 1) + 1
                % Here, we handle edge case of having one more datapoint
                % We just ignore the final datapoint
                output_struct.v(:, :, v_idx) = dataSeg(1:(end-1), :);
            else
                % Otherwise, we are missing too much data and need to store NaNs
                output_struct.v(:, :, v_idx) = ...
                    nan*ones(size(output_struct.v(:, :, v_idx)));
            end
        end
    end
end

% ---------- Take care of Neutral ------------- %
% Store neutral start and end times for segmentation
neutr_start = table2array(protTable(:, 'neutral_start')) - t_befStart;
neutr_end = neutr_start + t_befStart + t_neutr_seg + t_aftSeg;

for n_idx = 1:numN
    % If this participant had missing timings, then that means they ended
    % the protocol early; just store NaNs
    if n_idx > length(neutr_start)
        output_struct.n(:, :, n_idx) = ...
            nan*ones(size(output_struct.n(:, :, n_idx)));
    else
        % If a NaN is stored, ignore this timing and store NaNs
        if isnan(neutr_start(n_idx))
            output_struct.n(:, :, n_idx) = ...
                nan*ones(size(output_struct.n(:, :, n_idx)));
        else
            % Isolate the segment of data corresponding to this timing
            dataSeg = yArr(tVec >= neutr_start(n_idx) & ...
                tVec <= neutr_end(n_idx), :);
            
            % If we need a separate baseline for each response
            if sepBaselines
                % Compute baseline for this segment of data
                baseline = mean(yArr(...
                    tVec >= neutr_start(n_idx) + t_befStart - t_sepBaseline ...
                    & tVec <= neutr_start(n_idx) + t_befStart, :), 1);
                
                % Subtract away the baseline
                for feat_idx = 1:size(yArr, 2)
                    dataSeg(:, feat_idx) = dataSeg(:, feat_idx) - ...
                        baseline(feat_idx);
                end
            end
            
            % Make sure enough elements exist in this segment
            if size(dataSeg, 1) == size(output_struct.n, 1)
                % Store away
                output_struct.n(:, :, n_idx) = dataSeg;
            elseif size(dataSeg, 1) == size(output_struct.n, 1) - 1
                % Here, we handle edge case of having one less datapoint
                % We just repeat the final datapoint
                output_struct.n(:, :, n_idx) = [dataSeg; dataSeg(end, :)];
            elseif size(dataSeg, 1) == size(output_struct.n, 1) + 1
                % Here, we handle edge case of having one more datapoint
                % We just ignore the final datapoint
                output_struct.n(:, :, n_idx) = dataSeg(1:(end-1), :);
            else
                % Otherwise, we are missing too much data and need to store NaNs
                output_struct.n(:, :, n_idx) = ...
                    nan*ones(size(output_struct.n(:, :, n_idx)));
            end
        end
    end
end


% ---------- Take care of Trauma ------------- %
% Store trauma start and end times for segmentation
traum_start = table2array(protTable(:, 'trauma_start')) - t_befStart;
traum_end = traum_start + t_befStart + t_traum_seg + t_aftSeg;

for t_idx = 1:numT
    % If this participant had missing timings, then that means they ended
    % the protocol early; just store NaNs
    if t_idx > length(traum_start)
        output_struct.t(:, :, t_idx) = ...
            nan*ones(size(output_struct.t(:, :, t_idx)));
    else
        % If a NaN is stored, ignore this timing and store NaNs
        if isnan(traum_start(t_idx))
            output_struct.t(:, :, t_idx) = ...
                nan*ones(size(output_struct.t(:, :, t_idx)));
        else
            % Isolate the segment of data corresponding to this timing
            dataSeg = yArr(tVec >= traum_start(t_idx) & ...
                tVec <= traum_end(t_idx), :);
            
            % If we need a separate baseline for each response
            if sepBaselines
                % Compute baseline for this segment of data
                baseline = mean(yArr(...
                    tVec >= traum_start(t_idx) + t_befStart - t_sepBaseline ...
                    & tVec <= traum_start(t_idx) + t_befStart, :), 1);
                
                % Subtract away the baseline
                for feat_idx = 1:size(yArr, 2)
                    dataSeg(:, feat_idx) = dataSeg(:, feat_idx) - ...
                        baseline(feat_idx);
                end
            end
            
            % Make sure enough elements exist in this segment
            if size(dataSeg, 1) == size(output_struct.t, 1)
                % Store away
                output_struct.t(:, :, t_idx) = dataSeg;
            elseif size(dataSeg, 1) == size(output_struct.t, 1) - 1
                % Here, we handle edge case of having one less datapoint
                % We just repeat the final datapoint
                output_struct.t(:, :, t_idx) = [dataSeg; dataSeg(end, :)];
            elseif size(dataSeg, 1) == size(output_struct.t, 1) + 1
                % Here, we handle edge case of having one more datapoint
                % We just ignore the final datapoint
                output_struct.t(:, :, t_idx) = dataSeg(1:(end-1), :);
            else
                % Otherwise, we are missing too much data and need to store NaNs
                output_struct.t(:, :, t_idx) = ...
                    nan*ones(size(output_struct.t(:, :, t_idx)));
            end
        end
    end
end



end

