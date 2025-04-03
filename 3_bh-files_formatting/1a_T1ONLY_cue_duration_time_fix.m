% Specify the folder where the .mat files are located
folder = '/data/project/STAR/STAR_preprocessed/L_OnsDur/'; 

% Get a list of all .mat files in the folder
matFiles = dir(fullfile(folder, '*.mat'));

% Loop through each .mat file
for k = 1:length(matFiles)
    % Get the file name
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    
    % Load the .mat file
    load(fullFileName);  % load variables from the .mat file into the workspace
    
    % Check if 'durations' exists in the loaded variables
    if exist('durations', 'var')
        % Modify cell 1 and cell 4 in 'durations'
        durations{1} = 2 * ones(size(durations{1}));
        durations{4} = 2 * ones(size(durations{4}));
        
        % Save the modified data back to the same .mat file
        save(fullFileName, 'names', 'onsets', 'durations' );  % Overwrite the old file
        
        % Clear 'durations' variable to avoid conflict with the next file
        clear durations;
    else
        % Display a warning if 'durations' is not found in the file
        warning('Variable "durations" not found in file: %s', fullFileName);
    end
end
