% filepath: /Users/z/Documents/2_STAR_MRI/scripts/mri_preprocessing/4_analysis/ROI_coverage.m
%-----------------------------------------------------------------------
% Voxel overlap analysis for ROIs
%-----------------------------------------------------------------------

clear; clc;

dataroot = '/data/project/STAR/STAR_preprocessed';  % Preprocessed data directory
subs = {  %'sub-ID01', 'sub-ID02', 'sub-ID03','sub-ID04', ... % London pts 
%   'sub-ID05', 'sub-ID06', 'sub-ID07', 'sub-ID08', ...
%  'sub-ID09', 'sub-ID10' , 'sub-ID11', 'sub-ID12',...
%   'sub-ID13', 'sub-ID14', 'sub-ID15', 'sub-ID16', ...
%   'sub-ID17','sub-ID18',... 
%   'sub-ID19','sub-ID20', ...
%   'sub-ID21','sub-ID22',... 
%   'sub-ID23','sub-ID24', ...
%   'sub-ID25', 'sub-ID26', 'sub-ID62','sub-ID63', ...
%   'sub-ID64', 'sub-ID65', 'sub-ID66', 'sub-ID67', ...
%   'sub-ID68', 'sub-ID69', 'sub-ID70', 'sub-ID71', ...
%   'sub-ID72','sub-ID73', ...
%   'sub-ID74'
% 'sub-ID47', 'sub-ID48', 'sub-ID49', 'sub-ID50', ... % Newcastle pts
% 'sub-ID51', 'sub-ID52', 'sub-ID53', 'sub-ID54', ...
% 'sub-ID55', 'sub-ID56', 'sub-ID57', 'sub-ID58', ...
% 'sub-ID59', 'sub-ID60', 'sub-ID61', 'sub-ID62' 
'sub-ID27', 'sub-ID28', 'sub-ID29', 'sub-ID30', ... % Manchester pts
'sub-ID31', 'sub-ID32', 'sub-ID33', 'sub-ID34', ...
'sub-ID35', 'sub-ID36', 'sub-ID37', 'sub-ID38', ...
'sub-ID39', 'sub-ID40', 'sub-ID41', 'sub-ID42', ...
'sub-ID43', 'sub-ID44', 'sub-ID45', 'sub-ID46'
};    % List of subject IDs

% Preallocate results cell array
results = {};

% Loop over subjects and compute mask voxel count for each
results = {};  % Preallocate results cell array

for i = 1:length(subs)
    subj = subs{i};
    fprintf('Processing subject: %s\n', subj);
    
    % Define subject mask path (using fullfile for portability)
    mask_path = fullfile(dataroot, 'M_fmriprep', 'data', 'derivatives', subj, 'ses-mri01', 'func', 'analysis', 'mask.nii');
    if exist(mask_path, 'file') ~= 2
        warning('Full brain mask not found for %s, skipping.', subj);
        continue;
    end
    
    % Load subject full brain mask using SPM functions
    subj_nii = spm_vol(mask_path);
    subj_data = spm_read_vols(subj_nii);
    
    % Create binary mask and count the voxels in the mask
    subj_bin = subj_data > 0;
    voxel_count = sum(subj_bin(:));
    
    % Append results to cell array: subject ID and voxel count
    results = [results; {subj, voxel_count}];
end

% Create table from results and save to CSV
T = cell2table(results, 'VariableNames', {'SubjectID', 'ROI', 'OverlapVoxels', 'ROIVoxels', 'ProportionOverlap'});
disp(T);
csv_file = fullfile(dataroot, 'roi_overlap_all_subjects_M1.csv');
writetable(T, csv_file);
fprintf('Results saved to %s\n', csv_file);