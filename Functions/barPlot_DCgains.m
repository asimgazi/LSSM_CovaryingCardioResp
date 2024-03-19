function barPlot_DCgains(plotTitle, fileName, aParams, Cparams, featNames, ...
    inputsOpt, varargin)
% This function is passed parameters from a C matrix (state-to-output)
% terms and terms from a B matrix (input-to-state terms) and terms from an
% A matrix (state dynamics) to get the DC gain C((I-A)^{-1})B
% and creates bar plots to compare the product's values between the various
% inputs

% !!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!
% Note: I didn't bother changing too many of the comments and names from the
% barPlot_CBproducts function, so many of the variable names and comments
% may not match up

% ---------- Inputs ---------- %
% plotTitle: Just so the plots have some meaningful title
% fileName: file name to use when saving figure
% aParams: either a cell array of A matrices, one per subject, or a single
%          A matrix for a single subject
% Cparams: either a cell array of C matrices, one per subject, or a single
%          C matrix for a single subject
% featNames: Cell array of feature names to use when labeling
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

% I already know that I am ignoring legend entries, so suppress this wrning
warning('off', 'MATLAB:legend:IgnoringExtraEntries')

% Check if inputs is a cell array or not; if cell array, that means we need
% to create plot for each element of the product across entire group
if iscell(Cparams)
    % To know how big to initialize my array to store everything, I need to
    % know the number of inputs and the number of features, which is the
    % first dimension of C
    
    % Initialize an array to store the all the input terms
    % I am going to transpose the products to have each input as its own
    % row to make things more easily mappable to bar plots with x axis as
    % elements
    % I am going to stack the subjects' inputs in 3rd dimension
    CBproducts = zeros(numInputs, size(Cparams{1}, 1), numel(Cparams));
    
    % Compute products and stack
    for sub = 1:numel(Cparams)
        CBproducts(:, :, sub) = ...
            transpose(Cparams{sub}*...
            inv(eye(size(aParams{sub}, 1)) - aParams{sub})*...
            inputs{sub});
    end
    
    % Create mean +- SEM figure
    fig = figure();
    set(fig,'Visible','on')
    % colStrings = {'m', 'r', 'g'};   % Colors
    colStrings = {'g', 'm', 'r'};   % Colors
    
    % Now, I need to compute means and SEMs across subjects,
    % element by element
    mean_CBproducts = mean(CBproducts, 3);
    std_CBproducts = std(CBproducts, 0, 3);
    % For SEM
    plusMinus_CBproducts = (1/sqrt(numel(Cparams)))*std_CBproducts;


    % Let's switch things around so I have neutral first, then stimulation,
    % then trauma
    if numInputs > 2
        dummy = [mean_CBproducts(3, :); mean_CBproducts(1, :); ...
            mean_CBproducts(2, :)];
        mean_CBproducts = dummy;

        dummy = [std_CBproducts(3, :); std_CBproducts(1, :); ...
            std_CBproducts(2, :)];
        std_CBproducts = dummy;

        dummy = [plusMinus_CBproducts(3, :); plusMinus_CBproducts(1, :); ...
            plusMinus_CBproducts(2, :)];
        plusMinus_CBproducts = dummy;
    end


    % Now, I can create overlaid plot
    baseX = 1:size(mean_CBproducts, 2);   % For x-axis
    hold on
    for inp = 1:numInputs
        % Horizontal shift is so data can be seen better
        errorbar(baseX + (inp-numInputs+1)*0.1, mean_CBproducts(inp, :), ...
            plusMinus_CBproducts(inp, :), [colStrings{inp}, 'o']);
    end
    hold off
    
    % Create legend
    if numInputs == 1
        legend({'DCgain_v'});
    elseif numInputs == 2
        if strcmp(inputsOpt, 'vn')
            legend({'DCgain_v', 'DCgain_n'});
        else
            legend({'DCgain_v', 'DCgain_t'});
        end
    else
        % legend({'DCgain_v', 'DCgain_t', 'DCgain_n'});
        legend({'DCgain_n', 'DCgain_s', 'DCgain_t'});
    end
    
    % Label x axis
    xticks(1:size(CBproducts, 2));
    xticklabels(featNames);
    xtickangle(45);
    xlabel('Feature')
    
    % Label y axis
    ylabel('Value')
    
    % Title the plot and save
    title([plotTitle, ' - Mean \pm SEM'])
    savefig(fig, [fileName, '_meanSEM']);
    close all
    
    % Create boxplot figure
    fig = figure();
    set(fig,'Visible','on')
    
    % I want to reorganize things to make the application of MATLAB's
    % built-in boxplot function easier to do
    
    % Initialize an array with the total possible elements
    boxplot_Ready = zeros(size(CBproducts, 2)*numel(inputs), 2, numInputs);
    
    % Loop through inputs first
    for inp = 1:numInputs
        % Initialize a row counter
        rowCount = 1;
        
        % Now, loop through elements
        for dim = 1:size(CBproducts, 2)
            % Loop through subjects
            for sub = 1:numel(inputs)
                % Store with corresponding element number (i.e.,
                % grouping variable)
                boxplot_Ready(rowCount, :, inp) = ...
                    [dim, CBproducts(inp, dim, sub)];
                rowCount = rowCount + 1;    % Increment
            end
        end
    end
    
    % Now I can create overlaid plot
    baseX = 1:size(CBproducts, 2);   % For x-axis positions
    outlierSymbs = {'m+', 'ro', 'gx'};  % Just to make outliers different
    
    % Loop through inputs
    hold on
    for inp = 1:numInputs
        % Horizontal shift is so data can be seen better
        boxplot(boxplot_Ready(:, 2, inp), boxplot_Ready(:, 1, inp), ...
            'Positions', baseX + (inp-numInputs+1)*0.1, ...
            'Colors', colStrings{inp}, 'Symbol', outlierSymbs{inp});
    end
    
    
    % Create legend, but leave empty for now
    legend('')
    % Add extra lines to make legending easy
    if numInputs == 1
        hold on
        plot([nan, nan], [nan, nan], colStrings{1}, 'DisplayName', 'DCgain_v');
        hold off
    elseif numInputs == 2
        if strcmp(inputsOpt, 'vn')
            plot([nan, nan], [nan, nan], colStrings{1}, 'DisplayName', 'DCgain_v');
            plot([nan, nan], [nan, nan], colStrings{2}, 'DisplayName', 'DCgain_n');
        else
            plot([nan, nan], [nan, nan], colStrings{1}, 'DisplayName', 'DCgain_v');
            plot([nan, nan], [nan, nan], colStrings{2}, 'DisplayName', 'DCgain_t');
        end
    else
        % plot([nan, nan], [nan, nan], colStrings{1}, 'DisplayName', 'DCgain_v');
        plot([nan, nan], [nan, nan], colStrings{1}, 'DisplayName', 'DCgain_s');
        plot([nan, nan], [nan, nan], colStrings{2}, 'DisplayName', 'DCgain_t');
        plot([nan, nan], [nan, nan], colStrings{3}, 'DisplayName', 'DCgain_n');
    end
    hold off
    
    % Label x axis
    xticks(1:size(CBproducts, 2));
    xticklabels(featNames);
    xtickangle(45);
    xlabel('Feature')
    
    % Label y axis
    ylabel('Value')
    
    % Title the plot and save
    title([plotTitle, '- Box and Whiskers'])
    savefig(fig, [fileName, '_boxWhiskers']);
    
    % Close all figures
    close all
else
    % Create one figure
    fig = figure();
    set(fig,'Visible','on')
    
    % I am going to transpose the product to have each input's set of terms
    % stored in a separate row 
    CBproduct = transpose(Cparams*inv(eye(size(aParams, 1)) - aParams)*inputs);
    
    % Now, I can create overlaid barplot, where widths are halved each time
    hold on
    for inp = 1:numInputs
        bar(CBproduct(inp, :), 0.5^(inp-1))
    end
    hold off
    
    % Create legend
    if numInputs == 1
        legend({'DCgain_v'});
    elseif numInputs == 2
        if strcmp(inputsOpt, 'vn')
            legend({'DCgain_v', 'DCgain_n'});
        else
            legend({'DCgain_v', 'DCgain_t'});
        end
    else
        legend({'DCgain_v', 'DCgain_t', 'DCgain_n'});
    end
    
    % Label x axis
    xticks(1:size(CBproduct, 2));
    xticklabels(featNames);
    xtickangle(45);
    xlabel('Feature')
    
    % Label y axis
    ylabel('Value')
    
    % Title the plot, save, and close
    title(plotTitle)
    savefig(fig, fileName);
    close all
end

end

