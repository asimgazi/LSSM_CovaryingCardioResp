function norm_CBproducts(plotTitle, fileName, Cparams, ...
    inputsOpt, varargin)
% This function is passed parameters from a C matrix (state-to-output)
% terms and terms from a B matrix (input-to-state terms)
% and compares the norms of the input-specific C*b products

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
% ---> the next varargin input needs to be the b parameters arranged into
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


% Create one figure
fig = figure();
set(fig,'Visible','on')

% Check if inputs is a cell array or not; if cell array, that means we need
% to create plot across entire group
if iscell(Cparams)
    % Define a color array to try to have different colors for all of the
    % different subjects' specific datapoints
    % subjColors = (1/255)*[255, 85, 85; ...
    %     211, 95, 95; ...
    %     255, 153, 85; ...
    %     211, 141, 95; ...
    %     172, 157, 147; ...
    %     255, 221, 85; ...
    %     211, 188, 95; ...
    %     221, 255, 85; ...
    %     188, 211, 95; ...
    %     153, 255, 85; ...
    %     141, 211, 95; ...
    %     124, 145, 111; ...
    %     85, 255, 85; ...
    %     95, 211, 95; ...
    %     85, 255, 153; ...
    %     95, 211, 141; ...
    %     85, 255, 221; ...
    %     95, 211, 188; ...
    %     85, 221, 255; ...
    %     95, 188, 211; ...
    %     85, 153, 255; ...
    %     85, 85, 255; ...
    %     153, 85, 255; ...
    %     221, 85, 255; ...
    %     255, 85, 221; ...
    %     211, 95, 188; ...
    %     255, 85, 153];
    
    % We are just going to stick with light gray
    subjColors = (210/255)*ones(numel(Cparams), 3);
    
    if numel(Cparams) > size(subjColors, 1)
        % Fill rest randomly
        subjColors = [subjColors; ...
            rand(numel(Cparams) - size(subjColors, 1), 3)];
    end
    
    % Create a new array to store all the 2-norms
    norms_all = zeros(numel(Cparams), numInputs);
    
    % Loop through all and compute norms
    for sub = 1:numel(Cparams)
        % Compute norms for each C*b product
        for inp = 1:numInputs
            norms_all(sub, inp) = norm(Cparams{sub}*inputs{sub}(:, inp));
        end
    end


    % WE ARE GOING TO REORDER THE COLUMNS SUCH THAT IT IS NEUTRAL FIRST,
    % THEN STIMULATION, THEN TRAUMA
    if numInputs > 2
        dummy = [norms_all(:, 3), norms_all(:, 1), norms_all(:, 2)];
        norms_all = dummy;
    end
    
    % Compute p values for pairwise comparisons between groups, if applic.
    if numInputs == 2
        % Either signed rank test or paired t test
        if swtest(norms_all(:, 1)) || swtest(norms_all(:, 2))
            [pVal1, ~] = signrank(norms_all(:, 1), norms_all(:, 2));
            statsAppend = [': P_1 = ', num2str(pVal1), '(SR)'];
        else
            [~, pVal1] = ttest(norms_all(:, 1), norms_all(:, 2));
            statsAppend = [': P_1 = ', num2str(pVal1), '(TT)'];
        end
    elseif numInputs == 3
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
    boxplot(norms_all);
    
    hold on
    % Scatter points and connect the ones from the same subject with lines
    x_coord = 1:numInputs;
    for sub = 1:numel(inputs)
        p = plot(x_coord, norms_all(sub, :), '-o', ...
            'Color', subjColors(sub, :));
        p.MarkerFaceColor = subjColors(sub, :);
        p.MarkerSize = 4;
    end
    hold off
    
    % Label x axis
    xticks(1:numInputs);
    if numInputs == 1
        xticklabels({'C*b_v'});
    elseif numInputs == 2
        if strcmp(inputsOpt, 'vn')
            xticklabels({'C*b_v', 'C*b_n'});
        else
            xticklabels({'C*b_v', 'C*b_t'});
        end
    else
        % xticklabels({'C*b_v', 'C*b_t', 'C*b_n'});
        xticklabels({'C*b_n', 'C*b_s', 'C*b_t'});
    end
    
    % Label y axis
    ylabel('Magnitude')
    
    % Overall title, save, and close
    if numInputs > 1
        title([plotTitle, statsAppend])
    else
        title(plotTitle);
    end
    
    savefig(fig, fileName);
    close all
else
    % Initialize vector of norms
    norms = zeros(1, numInputs);
    
    % Compute norms
    for inp = 1:numInputs
        norms(inp) = norm(Cparams*inputs(:, inp));
    end
    
    % Create a bar plot to compare these norms
    bar(norms);
    
    % Label x axis
    xticks(1:numInputs);
    if numInputs == 1
        xticklabels({'C*b_v'});
    elseif numInputs == 2
        if strcmp(inputsOpt, 'vn')
            xticklabels({'C*b_v', 'C*b_n'});
        else
            xticklabels({'C*b_v', 'C*b_t'});
        end
    else
        xticklabels({'C*b_v', 'C*b_t', 'C*b_n'});
    end
    
    % Label y axis
    ylabel('Magnitude')
    
    % Title the plot, save, and close
    title(plotTitle);
    savefig(fig, fileName);
    close all
end

end

