function [timeTrain, timeTest, varargout] = split_trainTest(timeVec, ...
    protTimings, trTest_split, varargin)
% This function takes a time vector, protocol timings, and a splitting
% variable that dictates where to split the data into training and testing.
% Then, for any arrays provided to the function, it splits into train and
% test accordingly and returns them (along with time, of course)

% --------- Inputs ---------- %
% timeVec: Nx1 array of timepoints
% protTimings: table of start and end times for each protocol condition
% trTest_split: scalar value that dictates where to split data (e.g., if
%               trTest_split == 6.5, that means to split the data in
%               between the end of condition 6 and the beginning of
%               condition 7; note that trTest_split > 0 and < 11
% (Optional)
% Additional arrays that are input are assumed to be Nx___ arrays, where
% each array will be split in the same way that timeVec gets split into
% timeTrain and timeTest


% --------- Outputs --------- %
% timeTrain: 1-D array that stores the timepoints corresponding to training
%            data
% timeTest: 1-D array that stores the timepoints corresponding to testing
%           data
% (Optional)
% Additional pairs of arrays corresponding to training and testing splits
% of the optional arrays input to the function


% This function is created for the DARPA PTSD dataset, so the trTest_split
% convention follows the following:
% |N1|N2|T1V1|T2V2|V3|V4|N3|N4|T3V5|T4V6| (10 full conditions in total)
% So, if a user inputs trTest_split == 6.5, the split will occur after V4
% before N3
if trTest_split < 1
    startTime = timeVec(1);
    endTime = table2array(protTimings(1, 'neutral_start'));
elseif trTest_split > 1 && trTest_split < 2
    startTime = table2array(protTimings(1, 'neutral_end'));
    endTime = table2array(protTimings(2, 'neutral_start'));
elseif trTest_split > 2 && trTest_split < 3
    startTime = table2array(protTimings(2, 'neutral_end'));
    endTime = table2array(protTimings(1, 'trauma_start'));
elseif trTest_split > 3 && trTest_split < 4
    startTime = table2array(protTimings(1, 'vns_end'));
    endTime = table2array(protTimings(2, 'trauma_start'));
elseif trTest_split > 4 && trTest_split < 5
    startTime = table2array(protTimings(2, 'vns_end'));
    endTime = table2array(protTimings(3, 'vns_start'));
elseif trTest_split > 5 && trTest_split < 6
    startTime = table2array(protTimings(3, 'vns_end'));
    endTime = table2array(protTimings(4, 'vns_start'));
elseif trTest_split > 6 && trTest_split < 7
    startTime = table2array(protTimings(4, 'vns_end'));
    endTime = table2array(protTimings(3, 'neutral_start'));
elseif trTest_split > 7 && trTest_split < 8
    startTime = table2array(protTimings(3, 'neutral_end'));
    endTime = table2array(protTimings(4, 'neutral_start'));
elseif trTest_split > 8 && trTest_split < 9
    startTime = table2array(protTimings(4, 'neutral_end'));
    endTime = table2array(protTimings(3, 'trauma_start'));
elseif trTest_split > 9 && trTest_split < 10
    startTime = table2array(protTimings(5, 'vns_end'));
    endTime = table2array(protTimings(4, 'trauma_start'));
else
    startTime = table2array(protTimings(6, 'vns_end'));
    endTime = timeVec(end);
end

% Use the fractional component of trTest_split to figure out exactly where
% in between does the split need to occur
frComp = trTest_split - floor(trTest_split);
splitTime = frComp*(endTime - startTime) + startTime;

% Split timeVec first
timeTrain = timeVec(timeVec < splitTime);
timeTest = timeVec(timeVec >= splitTime);

% Initialize varargout, knowing there will be twice as many arrays
varargout = cell(1, 2*numel(varargin));

% Now, split the varargin into train and test
for in = 1:numel(varargin)
    % Train
    varargout{2*in-1} = varargin{in}(timeVec < splitTime, :);
    
    % Test
    varargout{2*in} = varargin{in}(timeVec >= splitTime, :);
end

end

