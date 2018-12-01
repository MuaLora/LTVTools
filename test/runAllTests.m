function results = runAllTests()
%% Run all tests

% Automatically find all the folders and call 'runall.m' for each
allFiles = dir(testroot);
cellArray = {allFiles.name};
allFolders = cellArray([allFiles.isdir]);
folderNames = allFolders(3:end);
results = cell(numel(folderNames),1);
passedFlag = true;

% Run all the folders
for i = 1:numel(folderNames)
    fprintf('=========================================================\n');
    fprintf('### Running Test from the following directory\n%s\n',...
        fullfile(testroot,folderNames{i}));
    fprintf('=========================================================\n');
    results{i} = fevalin(fullfile(testroot,folderNames{i}),'runall');
    
    % If all the tests did not pass then raise the flag
    if ~all([results{i}.Passed])
        passedFlag = false;
    end
end

% Print Status if there are failures or not
if ~passedFlag
    fprintf('### There were failures in the test results.\n');
else
    fprintf('### There were no test failures.\n');
end
end

%% Evaluate function at specific location
function out = fevalin(location, fcn)
oldfolder = cd(location);
out = feval(fcn);
cd(oldfolder);
end