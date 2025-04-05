
dataroot = '/data/project/STAR/STAR_preprocessed'; % Preprocessed data directory
roiroot = '/data/project/STAR/STAR_preprocessed/ROIanalysis_20250106/';
% Define ROI mask filenames
roi_files = {'roi1.nii', 'roi2.nii', 'roi3.nii'};

% Define subject IDs (update as needed)
subject_nums = 1:2;  % e.g., sub-ID01 and sub-ID02
n_subjects = length(subject_nums);

% Preallocate results
results = [];

% Loop over subjects
for s = 1:n_subjects
    % Format subject ID with leading zeros
    subj_id = sprintf('ID%02d', subject_nums(s));
    mask_path = fullfile('/data/project/STAR/STAR_preprocessed/L_fmriprep/data/derivatives', ...
        ['sub-' subj_id], 'ses-mri01', 'func', 'mask.nii');

    % Load subject's full brain mask
    if exist(mask_path, 'file') ~= 2
        warning('Mask not found for %s, skipping.', subj_id);
        continue;
    end

    subj_nii = spm_vol(mask_path);
    subj_data = spm_read_vols(subj_nii);
    subj_bin = subj_data > 0;

    % Loop over ROIs
    for r = 1:length(roi_files)
        roi_path = roi_files{r};

        if exist(roi_path, 'file') ~= 2
            warning('ROI file %s not found, skipping.', roi_path);
            continue;
        end

        % Load and binarize ROI
        roi_nii = spm_vol(roi_path);
        roi_data = spm_read_vols(roi_nii);
        roi_bin = roi_data > 0;

        % Calculate overlap
        overlap = roi_bin & subj_bin;
        overlap_voxels = sum(overlap(:));
        roi_voxels = sum(roi_bin(:));
        prop_overlap = overlap_voxels / roi_voxels;

        % Store results
        results = [results; 
            {['sub-' subj_id], roi_files{r}, overlap_voxels, roi_voxels, prop_overlap}];
    end
end

% Create table
T = cell2table(results, ...
    'VariableNames', {'SubjectID', 'ROI', 'OverlapVoxels', 'ROIVoxels', 'ProportionOverlap'});

% Display and save
disp(T);
writetable(T, 'roi_overlap_all_subjects.csv');
