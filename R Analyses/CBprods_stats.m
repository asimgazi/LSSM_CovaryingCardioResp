% Code for CB products to compute statistics / prep for R 

% It is assumed you have a variable that stores all the CB products you are
% interested in, with first dimension being number of inputs, second
% dimenion being number of features, and third dimension being number of
% subjects
CBprods = CBproducts_all;

%% First, we will compute shapiro-wilk normality statistics
% We will create a table we can easily copy and paste into excel
numFeats = size(CBprods, 2);
numInputs = size(CBprods, 1);
numSubj = size(CBprods, 3);

% Initialize table
swtest_results = -1*ones(numInputs, numFeats);

% Loop through and fill table
for feat = 1:numFeats
    for input = 1:numInputs
        dataVec = CBprods(input, feat, :);
        swtest_results(input, feat) = swtest(dataVec);
    end
end

%% Now, we will create a CSV file that we can load into R for further anlalysis
% Since R is more reliable for statistics, I will just resort to outputting
% a CSV file and then do the rest of the statistical comparisons in R

% Looks like matrix format where input is column (the factor that I am
% comparing) and rows can be subjects. What I can do is just create a
% separate CSV for each feature to make things easier to load and code in R

% I will reorganize CB products to have features be the third dimension,
% inputs the second, and subjects the first
CBprods_csvReady = zeros(numSubj, numInputs, numFeats);

% Reorganize and store
for feat = 1:numFeats
    for subj = 1:numSubj
        for input = 1:numInputs
            CBprods_csvReady(subj, input, feat) = CBprods(input, feat, subj);
        end
    end
end

% Loop through and use feature names stored in featNames_desire to label
% CSV files
for feat = 1:numFeats
    writematrix(CBprods_csvReady(:, :, feat), ...
        ['CBprods_csvReady_', featNames_desir{feat}, '.csv']);
end
