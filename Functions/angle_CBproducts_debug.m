function angle_CBproducts_debug(plotTitle, fileName, Cparams, inputsOpt, varargin)
% This function is passed parameters from a C matrix (state-to-output)
% terms and terms from a B matrix (input-to-state terms)
% and compares the angles between input-specific terms

% ---------- Inputs ---------- %
% plotTitle: Just so the plots have some meaningful title
% fileName: file name to use when saving figure
% Cparams: either a cell array of C matrices, one per subject, or a single
%          C matrix for a single subject
% inputsOpt: string that specifies which protocol conditions
%            inputs were created for - 
%            'v': just VNS
%            'vt': VNS and trauma recall
%            'vn': VNS and neutral condition
%            'vtn': VNS, trauma recall, and neutral
% ---> the next varargin inputs needs to be the b parameters arranged into
%      columns according to inputsOpt (e.g., if inputsOpt is 'vn', then the
%      first column needs to be b_v, and the 2nd column needs to be b_n)
% (optional) 
% 'elementNums' - if the user enters this flag, the next varargin needs to
%                 be a vector of element numbers to include in the analysis
%                 (useful to compare angles in lower dimension of one's
%                  choosing)

% Parse varargin
if strcmp(inputsOpt, 'vtn')
    numInputs = 3;  % Store number of inputs
    inputs = varargin{1};
elseif strcmp(inputsOpt, 'vt') || strcmp(inputsOpt, 'vn')
    numInputs = 2;  % Store number of inputs
    inputs = varargin{1};
else
    numInputs = 1;  % Store number of inputs
    inputs = varargin{1};
end
for arg = 2:length(varargin)
    if strcmp(varargin{arg}, 'elementNums')
        elementNums = varargin{arg + 1};
    end
end

% Set defaults
if ~exist('elementNums', 'var')
    if iscell(inputs)
        % Use all elements
        elementNums = 1:size(Cparams{1}*inputs{1}, 1);
    else
        % Use all elements
        elementNums = 1:size(Cparams*inputs, 1);
    end
end

% Compute number of pairs
if numInputs == 2
    numPairs = 1;
elseif numInputs == 3
    numPairs = 3;
else
    numPairs = 0;
end

% Check if inputs is a cell array or not; if cell array, that means we need
% to create a boxplot to plot an entire group's data
if iscell(inputs)
    % Create a new array to store all the angles
    angles_all = zeros(numel(inputs), numPairs);
    
    % Loop through all and compute angles, if there are pairs to work with
    if ~isempty(angles_all)
        for sub = 1:numel(inputs)
            % Compute angles between normalized pairs
            p = 1;  % Pair index
            for p1 = 1:(numInputs-1)
                for p2 = (p1+1):numInputs
                    % Compute products and store away desired elements
                    vec1_full = Cparams{sub}*inputs{sub}(:, p1);
                    vec1_desired = vec1_full(elementNums);
                    
                    vec2_full = Cparams{sub}*inputs{sub}(:, p2);
                    vec2_desired = vec2_full(elementNums);
                    
                    % Make vectors unit norm
                    vec1_normed = (1/norm(vec1_desired))*vec1_desired;
                    vec2_normed = (1/norm(vec2_desired))*vec2_desired;
                    
                    % Take arc cosine of the dot product and convert to degrees
                    angles_all(sub, p) = (180/pi)*acos(vec1_normed'*...
                        vec2_normed);
                    
                    p = p + 1;  % increment pair counter
                end
            end
        end
        
        % Compute p values for pairwise comparisons between groups, if applic.
        if numPairs == 2
            % Either signed rank test or paired t test
            if swtest(angles_all(:, 1)) || swtest(angles_all(:, 2))
                [pVal1, ~] = signrank(angles_all(:, 1), angles_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(angles_all(:, 1), angles_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(TT)'];
            end
        elseif numPairs == 3
            % First p value
            if swtest(angles_all(:, 1)) || swtest(angles_all(:, 2))
                [pVal1, ~] = signrank(angles_all(:, 1), angles_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(angles_all(:, 1), angles_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(TT)'];
            end
            % Second p value
            if swtest(angles_all(:, 3)) || swtest(angles_all(:, 2))
                [pVal2, ~] = signrank(angles_all(:, 3), angles_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(SR)'];
            else
                [~, pVal2] = ttest(angles_all(:, 3), angles_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(TT)'];
            end
            % Third p value
            if swtest(angles_all(:, 3)) || swtest(angles_all(:, 1))
                [pVal3, ~] = signrank(angles_all(:, 3), angles_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(SR)'];
            else
                [~, pVal3] = ttest(angles_all(:, 3), angles_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(TT)'];
            end
            
            statsAppend = [statsAppend_1, statsAppend_2, statsAppend_3];
        end
        
        % Create a boxplot for all norms, separated by input
        boxplot(angles_all);
        
        hold on
        % Scatter points and connect the ones from the same subject with lines
        x_coord = 1:numPairs;
        for sub = 1:numel(inputs)
            plot(x_coord, angles_all(sub, :), '-o')
        end
        hold off
        
        % Label x axis
        xticks(1:numPairs);
         % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'\theta_{C*b_v, C*b_n}'});
            else
                xticklabels({'\theta_{C*b_v, C*b_t}'});
            end
        else
            xticklabels({'\theta_{C*b_v, C*b_t}', '\theta_{C*b_v, C*b_n}', ...
                '\theta_{C*b_t, C*b_n}'});
        end
        
        % Label y axis
        ylabel('Angle (degrees)')
        
        % Overall title, save, and close
        if numPairs > 1
            title([plotTitle, statsAppend])
        else
            title(plotTitle);
        end
        
        % Save
        savefig(fig, fileName);
    end
    % Close
    close all
else
    angles = zeros(1, numPairs);
    
    % Compute angles between normalized pairs
    p = 1;  % Pair index
    for p1 = 1:(numInputs-1)
        for p2 = (p1+1):numInputs
            % Compute products and store away desired elements
            vec1_full = Cparams*inputs(:, p1);
            vec1_desired = vec1_full(elementNums);
            
            vec2_full = Cparams*inputs(:, p2);
            vec2_desired = vec2_full(elementNums);
            
            % Make vectors unit norm
            vec1_normed = (1/norm(vec1_desired))*vec1_desired;
            vec2_normed = (1/norm(vec2_desired))*vec2_desired;
            
            % Take arc cosine of the dot product and convert to degrees
            angles(p) = (180/pi)*acos(vec1_normed'*...
                vec2_normed);
            
            p = p + 1;  % increment pair counter
        end
    end
    
    % Create one figure
    fig = figure();
    set(fig,'Visible','on')
    
    subplot(1, 2, 1)
    % Create a bar plot to compare these angles
    bar(angles);
    
    % Label x axis
    xticks(1:numPairs);
    if numInputs == 2
        if strcmp(inputsOpt, 'vn')
            xticklabels({'\theta_{C*b_v, C*b_n}'});
        else
            xticklabels({'\theta_{C*b_v, C*b_t}'});
        end
    else
        xticklabels({'\theta_{C*b_v, C*b_t}', '\theta_{C*b_v, C*b_n}', ...
            '\theta_{C*b_t, C*b_n}'});
    end
    
    % Label y axis
    ylabel('Angle (degrees)')
    
    % Title the plot, save, and close
    title(plotTitle);
    
    subplot(1, 2, 2)
    % create a plot to display the three vectors
    hold on
    for inp = 1:numInputs
        vec_full = Cparams*inputs(:, inp);
        vec_desired = vec_full(elementNums);
        Xcoords = [0; vec_desired(1)];
        Ycoords = [0; vec_desired(2)];
        plot(Xcoords, Ycoords);
    end
    hold off
    xlabel('HR')
    ylabel('RR')
    legend({'C*b_v', 'C*b_t', 'C*b_n'})
    
    % Title the plot, save, and close
    title(plotTitle);
    
    savefig(fig, fileName);
    close all
end


end

