function normDiffs_CBproducts(plotTitle, fileName, Cparams, inputsOpt, ...
    varargin)
% This function is passed parameters from a C matrix (state-to-output)
% terms and terms from a B matrix (input-to-state terms)
% and compares the norms of the differences between input-specific 
% products

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
    % Create a new array to store all the norms
    norms_all = zeros(numel(inputs), numPairs);
    norms_div_all = zeros(numel(inputs), numPairs);
    norms_unit_all = zeros(numel(inputs), numPairs);
    
    % Loop through all and compute norms of diffs, if there are pairs
    if ~isempty(norms_all)
        for sub = 1:numel(inputs)
            p = 1;  % Pair index
            for p1 = 1:(numInputs-1)
                for p2 = (p1+1):numInputs
                    % Compute the norm of the pairwise difference
                    norms_all(sub, p) = norm(Cparams{sub}*inputs{sub}(:, p1) - ...
                        Cparams{sub}*inputs{sub}(:, p2));
                    norms_div_all(sub, p) = ...
                        norm(Cparams{sub}*inputs{sub}(:, p1) - ...
                        Cparams{sub}*inputs{sub}(:, p2))/...
                        (norm(Cparams{sub}*inputs{sub}(:, p1))*...
                        norm(Cparams{sub}*inputs{sub}(:, p2)));
                    norms_unit_all(sub, p) = norm(...
                        (Cparams{sub}*inputs{sub}(:, p1)/...
                        norm(Cparams{sub}*inputs{sub}(:, p1))) - ...
                        (Cparams{sub}*inputs{sub}(:, p2)/...
                        norm(Cparams{sub}*inputs{sub}(:, p2)))...
                        );
                    
                    p = p + 1;  % increment pair counter
                end
            end
        end
        
        % Compute p values for pairwise comparisons between groups, if applic.
        if numPairs == 2
            % Either signed rank test or paired t test
            if swtest(norms_all(:, 1)) || swtest(norms_all(:, 2))
                [pVal1, ~] = signrank(norms_all(:, 1), norms_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_all(:, 1), norms_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(TT)'];
            end
        elseif numPairs == 3
            % First p value
            if swtest(norms_all(:, 1)) || swtest(norms_all(:, 2))
                [pVal1, ~] = signrank(norms_all(:, 1), norms_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_all(:, 1), norms_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(TT)'];
            end
            % Second p value
            if swtest(norms_all(:, 3)) || swtest(norms_all(:, 2))
                [pVal2, ~] = signrank(norms_all(:, 3), norms_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(SR)'];
            else
                [~, pVal2] = ttest(norms_all(:, 3), norms_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(TT)'];
            end
            % Third p value
            if swtest(norms_all(:, 3)) || swtest(norms_all(:, 1))
                [pVal3, ~] = signrank(norms_all(:, 3), norms_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(SR)'];
            else
                [~, pVal3] = ttest(norms_all(:, 3), norms_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(TT)'];
            end
            
            statsAppend = [statsAppend_1, statsAppend_2, statsAppend_3];
        end
        
        % Create a boxplot for all norms, separated by input
        % Create one figure
        fig = figure();
        set(fig,'Visible','on')
        boxplot(norms_all);
        
        hold on
        % Scatter points and connect the ones from the same subject with lines
        x_coord = 1:numPairs;
        for sub = 1:numel(inputs)
            plot(x_coord, norms_all(sub, :), '-o')
        end
        hold off
        
        % Label x axis
        xticks(1:numPairs);
         % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'C*b_v - C*b_n'});
            else
                xticklabels({'C*b_v - C*b_t'});
            end
        else
            xticklabels({'C*b_v - C*b_t', 'C*b_v - C*b_n', ...
                'C*b_t - C*b_n'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Overall title, save, and close
        if numPairs > 1
            title([plotTitle, statsAppend])
        else
            title(plotTitle);
        end
        
        % Save
        savefig(fig, fileName);
        % Close
        close all
        
        
        
        % ------- Replicate for norms_div ----------- %
        % Compute p values for pairwise comparisons between groups, if applic.
        if numPairs == 2
            % Either signed rank test or paired t test
            if swtest(norms_div_all(:, 1)) || swtest(norms_div_all(:, 2))
                [pVal1, ~] = signrank(norms_div_all(:, 1), norms_div_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_div_all(:, 1), norms_div_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(TT)'];
            end
        elseif numPairs == 3
            % First p value
            if swtest(norms_div_all(:, 1)) || swtest(norms_div_all(:, 2))
                [pVal1, ~] = signrank(norms_div_all(:, 1), norms_div_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_div_all(:, 1), norms_div_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(TT)'];
            end
            % Second p value
            if swtest(norms_div_all(:, 3)) || swtest(norms_div_all(:, 2))
                [pVal2, ~] = signrank(norms_div_all(:, 3), norms_div_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(SR)'];
            else
                [~, pVal2] = ttest(norms_div_all(:, 3), norms_div_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(TT)'];
            end
            % Third p value
            if swtest(norms_div_all(:, 3)) || swtest(norms_div_all(:, 1))
                [pVal3, ~] = signrank(norms_div_all(:, 3), norms_div_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(SR)'];
            else
                [~, pVal3] = ttest(norms_div_all(:, 3), norms_div_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(TT)'];
            end
            
            statsAppend = [statsAppend_1, statsAppend_2, statsAppend_3];
        end
        
        % Create a boxplot for all norms, separated by input
        % Create one figure
        fig = figure();
        set(fig,'Visible','on')
        boxplot(norms_div_all);
        
        hold on
        % Scatter points and connect the ones from the same subject with lines
        x_coord = 1:numPairs;
        for sub = 1:numel(inputs)
            plot(x_coord, norms_div_all(sub, :), '-o')
        end
        hold off
        
        % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'(C*b_v - C*b_n)/norms of each'});
            else
                xticklabels({'(C*b_v - C*b_t)/norms of each'});
            end
        else
            xticklabels({'(C*b_v - C*b_t)/norms of each', ...
                '(C*b_v - C*b_n)/norms of each', ...
                '(C*b_t - C*b_n)/norms of each'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Overall title, save, and close
        if numPairs > 1
            title([plotTitle, statsAppend])
        else
            title(plotTitle);
        end
        
        % Save
        savefig(fig, [fileName, '_divNorms']);
        % Close
        close all
        
        
        % ------- Replicate for norms_unit ----------- %
        % Compute p values for pairwise comparisons between groups, if applic.
        if numPairs == 2
            % Either signed rank test or paired t test
            if swtest(norms_unit_all(:, 1)) || swtest(norms_unit_all(:, 2))
                [pVal1, ~] = signrank(norms_unit_all(:, 1), norms_unit_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_unit_all(:, 1), norms_unit_all(:, 2));
                statsAppend = [': P_1 = ', num2str(pVal1), '(TT)'];
            end
        elseif numPairs == 3
            % First p value
            if swtest(norms_unit_all(:, 1)) || swtest(norms_unit_all(:, 2))
                [pVal1, ~] = signrank(norms_unit_all(:, 1), norms_unit_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(SR)'];
            else
                [~, pVal1] = ttest(norms_unit_all(:, 1), norms_unit_all(:, 2));
                statsAppend_1 = [': P_{12} = ', num2str(pVal1), '(TT)'];
            end
            % Second p value
            if swtest(norms_unit_all(:, 3)) || swtest(norms_unit_all(:, 2))
                [pVal2, ~] = signrank(norms_unit_all(:, 3), norms_unit_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(SR)'];
            else
                [~, pVal2] = ttest(norms_unit_all(:, 3), norms_unit_all(:, 2));
                statsAppend_2 = ['; P_{23} = ', num2str(pVal2), '(TT)'];
            end
            % Third p value
            if swtest(norms_unit_all(:, 3)) || swtest(norms_unit_all(:, 1))
                [pVal3, ~] = signrank(norms_unit_all(:, 3), norms_unit_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(SR)'];
            else
                [~, pVal3] = ttest(norms_unit_all(:, 3), norms_unit_all(:, 1));
                statsAppend_3 = ['; P_{31} = ', num2str(pVal3), '(TT)'];
            end
            
            statsAppend = [statsAppend_1, statsAppend_2, statsAppend_3];
        end
        
        % Create a boxplot for all norms, separated by input
        % Create one figure
        fig = figure();
        set(fig,'Visible','on')
        boxplot(norms_unit_all);
        
        hold on
        % Scatter points and connect the ones from the same subject with lines
        x_coord = 1:numPairs;
        for sub = 1:numel(inputs)
            plot(x_coord, norms_unit_all(sub, :), '-o')
        end
        hold off
        
        % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'unit(C*b_v) - unit(C*b_n)'});
            else
                xticklabels({'unit(C*b_v) - unit(C*b_t)'});
            end
        else
            xticklabels({'unit(C*b_v) - unit(C*b_t)', ...
                'unit(C*b_v) - unit(C*b_n)', ...
                'unit(C*b_t) - unit(C*b_n)'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Overall title, save, and close
        if numPairs > 1
            title([plotTitle, statsAppend])
        else
            title(plotTitle);
        end
        
        % Save
        savefig(fig, [fileName, '_diffUnits']);
        % Close
        close all
    end
else
    % Initialize vector of norms
    % (initialize additional arrays for exploration purposes)
    norms = zeros(1, numPairs);
    norms_div = zeros(1, numPairs);
    norms_unit = zeros(1, numPairs);
    
    % Make sure there are pairs to work with
    if ~isempty(norms)
        % Compute norms of pairwise differences
        p = 1;  % Pair index
        for p1 = 1:(numInputs-1)
            for p2 = (p1+1):numInputs
                % Compute the norm of the difference
                norms(p) = norm(Cparams*inputs(:, p1) - Cparams*inputs(:, p2));
                norms_div(p) = ...
                    norm(Cparams*inputs(:, p1) - Cparams*inputs(:, p2))/...
                    (norm(Cparams*inputs(:, p1))*norm(Cparams*inputs(:, p2)));
                norms_unit(p) = norm(...
                    (Cparams*inputs(:, p1)/norm(Cparams*inputs(:, p1))) - ...
                    (Cparams*inputs(:, p2)/norm(Cparams*inputs(:, p2)))...
                    );
                p = p + 1;  % increment pair counter
            end
        end
        
        % Create a bar plot to compare these angles
        fig = figure();
        set(fig,'Visible','on')
        bar(norms);
        
        % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'C*b_v - C*b_n'});
            else
                xticklabels({'C*b_v - C*b_t'});
            end
        else
            xticklabels({'C*b_v - C*b_t', 'C*b_v - C*b_n', ...
                'C*b_t - C*b_n'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Title the plot, save, and close
        title(plotTitle);
        savefig(fig, fileName);
        close all
        
        
        % ------------ Replicate for norms_div ------------ %
        % Create a bar plot to compare these angles
        fig = figure();
        set(fig,'Visible','on')
        bar(norms_div);
        
        % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'(C*b_v - C*b_n)/norms of each'});
            else
                xticklabels({'(C*b_v - C*b_t)/norms of each'});
            end
        else
            xticklabels({'(C*b_v - C*b_t)/norms of each', ...
                '(C*b_v - C*b_n)/norms of each', ...
                '(C*b_t - C*b_n)/norms of each'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Title the plot, save, and close
        title(plotTitle);
        savefig(fig, [fileName, '_divNorms']);
        close all
        
        
        % ------------ Replicate for norms_unit ------------ %
        % Create a bar plot to compare these angles
        fig = figure();
        set(fig,'Visible','on')
        bar(norms_unit);
        
        % Label x axis
        xticks(1:numPairs);
        if numInputs == 2
            if strcmp(inputsOpt, 'vn')
                xticklabels({'unit(C*b_v) - unit(C*b_n)'});
            else
                xticklabels({'unit(C*b_v) - unit(C*b_t)'});
            end
        else
            xticklabels({'unit(C*b_v) - unit(C*b_t)', ...
                'unit(C*b_v) - unit(C*b_n)', ...
                'unit(C*b_t) - unit(C*b_n)'});
        end
        
        % Label y axis
        ylabel('Magnitude')
        
        % Title the plot, save, and close
        title(plotTitle);
        savefig(fig, [fileName, '_diffUnits']);
        close all
    end
end


end

