function inputArray = ...
    create_inputs(timeVec, protocol_timings, VNS_binary, VNSamps, inputsOpt)
% This function takes a time vector, protocol information, and the desired
% options and creates an array of input(s), where each column is an input
% and each row is a timepoint

% --------- Inputs ------------ %
% timeVec: Nx1 array of time points at regularly sampled intervals
% protocol_timings: table of protocol timings corresponding to DARPA PTSD
%                   VNS dataset
% VNS_binary: boolean dictating whether to treat VNS as binary or not
% VNSamps: 8x1 or empty array of VNS amplitudes, one for each VNS admin.
%          (can be left empty if VNS_binary is true)
% inputsOpt: string that specifies which protocol conditions to create
%            inputs for -
%            'v': just VNS (default)
%            'vt': VNS and trauma recall
%            'vn': VNS and the neutral condition
%            'vtn': VNS, trauma recall, and neutral

% --------- Outputs ----------- %
% inputArray: Nx__ (blank is dictated by inputsOpt) array of input values,
%             where each input is its own column


% Since VNS is always included as the first column of the input array, we
% take care of it first
VNSinput = zeros(length(timeVec), 1);

% Store VNS start times (we do not need to store end times since VNS always
% takes 120 s)
vns_start = table2array(protocol_timings(:, 'vns_start'));
vns_end = vns_start + 120;

% If we want VNSinput to either be 0 or 1
if VNS_binary
    % Loop through VNS administrations
    for vns_idx = 1:length(vns_start)
        % Make sure it's not NaN
        if ~isnan(vns_start(vns_idx))
            % Based on the timepoints within vns start and end, set the
            % corresponding VNSinput values to 1
            VNSinput(timeVec >= vns_start(vns_idx) & timeVec <= vns_end(vns_idx)) = ...
                ones(...
                length(timeVec(timeVec >= vns_start(vns_idx) & ...
                timeVec <= vns_end(vns_idx))), 1);
        end
    end
    
else
    % Otherwise, we have to use the VNSamps information to create our input
    % We are going to max-min scale our VNS amplitudes to keep the maximum
    % VNS input to be 1 (min is 0)
    VNSamps = (1/max(VNSamps))*VNSamps;
    
    % Loop through VNS administrations
    for vns_idx = 1:length(vns_start)
        % Make sure it's not NaN
        if ~isnan(vns_start(vns_idx))
            % Based on the timepoints within vns start and end, set the
            % corresponding VNSinput values to 1
            VNSinput(timeVec >= vns_start(vns_idx) & timeVec <= vns_end(vns_idx)) = ...
                VNSamps(vns_idx)*ones(...
                length(timeVec(timeVec >= vns_start(vns_idx) & ...
                timeVec <= vns_end(vns_idx))), 1);
        end
    end
end


% Now we create the trauma input, if desired
if strcmp(inputsOpt, 'vt') || strcmp(inputsOpt, 'vtn')
    % Initialize trauma input
    traumInput = zeros(length(timeVec), 1);
    
    % Store trauma start and end times
    traum_start = table2array(protocol_timings(:, 'trauma_start'));
    traum_end = table2array(protocol_timings(:, 'trauma_end'));
    
    % Loop through trauma recalls
    for t_idx = 1:length(traum_start)
        % Make sure it's not NaN
        if ~isnan(traum_start(t_idx))
            % Based on the timepoints within trauma start and end, set the
            % corresponding traumInput values to 1
            traumInput(timeVec >= traum_start(t_idx) & timeVec <= traum_end(t_idx)) = ...
                ones(...
                length(timeVec(timeVec >= traum_start(t_idx) & ...
                timeVec <= traum_end(t_idx))), 1);
        end
    end
end


% Now we create the neutral input, if desired
if strcmp(inputsOpt, 'vn') || strcmp(inputsOpt, 'vtn')
    % Initialize neutral input
    neutrInput = zeros(length(timeVec), 1);
    
    % Store trauma start and end times
    neutr_start = table2array(protocol_timings(:, 'neutral_start'));
    neutr_end = table2array(protocol_timings(:, 'neutral_end'));
    
    % Loop through trauma recalls
    for n_idx = 1:length(neutr_start)
        % Make sure it's not NaN
        if ~isnan(neutr_start(n_idx))
            % Based on the timepoints within neutral start and end, set the
            % corresponding neutrInput values to 1
            neutrInput(timeVec >= neutr_start(n_idx) & timeVec <= neutr_end(n_idx)) = ...
                ones(...
                length(timeVec(timeVec >= neutr_start(n_idx) & ...
                timeVec <= neutr_end(n_idx))), 1);
        end
    end
end


% Concatenate for final result
inputArray = VNSinput;

% Check if traumInput exists; if so, append
if exist('traumInput', 'var')
    inputArray = [inputArray, traumInput];
end

% Check if neutrInput exists; if so, append
if exist('neutrInput', 'var')
    inputArray = [inputArray, neutrInput];
end

end

