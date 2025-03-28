% Set the root folder
root = '/Users/z/Documents/MRI_data/scans_L';

% Initialize cell arrays to store data
onset = cell(1, numel(folders));
durations = cell(1, numel(folders));
trial_type = cell(1, numel(folders));
is_own = cell(1, numel(folders));
rating = cell(1, numel(folders));
response_time = cell(1, numel(folders));

% Get list of STAR participant numbers for parent folders
folders = dir(root);
folders = folders([folders.isdir]); % Filter out non-folders

% Loop through each participant folder
for i = 1:numel(folders)
    folder_path = fullfile(root, folders(i).name);

    % Load CSV files
    files = dir(fullfile(folder_path, 'T1/task', '*TraumaMemory*.csv'));
    for j = 1:numel(files)
        fullcsv = readtable(fullfile(files(j).folder, files(j).name));
        
        % Chop off unnecessary rows
        abbreviated = fullcsv(1:22, :);
        
        % Get onset times
        allonsets = table2array(abbreviated(:, {'TimeAtStartOfTrial', 'ElaborateCueOnset', 'PostElaborateFixOnset', 'RatingOnset', 'RestJitterOnset'}));
        onset = reshape(allonsets.', 1, []);
        
        % Make durations vector
        onset_times = onset(~isnan(onset)); % Remove NaNs
        durations = diff(onset_times);
        
        % Make trial_type vector
        trial_type = repmat({'cue', 'elaborate', 'pelaboratefix', 'rating', 'rest'}, 1, 16);
        trial_type = [trial_type, 'fixation'];
        
        % Make is_own and rating vectors
        is_own = repmat([0, 0, 0, 1, 0], 1, 16); 
        rating = nan(size(is_own));
        %response_time = nan(size(is_own));
       
        % Save data as .mat file
        save(sprintf('%s_run%d.mat', folders(i).name, j), 'onset', 'durations', 'trial_type', 'is_own', 'rating', 'response_time');
    end
end
