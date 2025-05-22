% filepath: /Users/z/Documents/2_STAR_MRI/scripts/mri_preprocessing/4_analysis/ROI_coverage.m
%-----------------------------------------------------------------------
% Voxel overlap analysis for ROIs
%-----------------------------------------------------------------------

clear; clc;

dataroot = '/data/project/STAR/STAR_preprocessed';  % Preprocessed data directory
roiroot = fullfile(dataroot, 'ROIanalysis_20250106'); % ROI analysis directory
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
roi_files = {'amygdala_wfu.nii', 'hippocampus_wfu.nii', 'vmpfc_bybox.nii'};   % ROI mask filenames
% once this works, change to have the left and right amygdala ROIs

% Preallocate results cell array
results = {};

% Loop over subjects
for i = 1:length(subs)
        subj = subs{i};
    fprintf('Processing subject: %s\n', subj);
    
    % Define subject mask path (using fullfile for portability)
    mask_path = fullfile(dataroot, 'M_fmriprep','data','derivatives', subj, 'ses-mri01', 'func', 'analysis', 'mask.nii'); %another layer down for mask for M
    if exist(mask_path, 'file') ~= 2
        warning('Full brain mask not found for %s, skipping.', subj);
        continue;
    end
    
    % Load subject full brain mask using SPM functions
    subj_nii = spm_vol(mask_path);
    subj_data = spm_read_vols(subj_nii);
    subj_bin = subj_data > 0;
    
    % Loop over each ROI
    for r = 1:length(roi_files)
        roi_path = fullfile(roiroot, roi_files{r});
        if exist(roi_path, 'file') ~= 2
            warning('ROI file %s not found for subject %s, skipping.', roi_files{r}, subj);
            continue;
        end
        
        % Load ROI mask
        roi_nii = spm_vol(roi_path);
        roi_data = spm_read_vols(roi_nii);
        roi_bin = roi_data > 0;
        
        % If dimensions do not match, reslice ROI to subject space
        if ~isequal(size(roi_bin), size(subj_bin))
            fprintf('Dimensions differ for ROI %s and subject %s. Reslicing ROI...\n', roi_files{r}, subj);
            % Create a two-row char matrix where first row is the reference and 
            % the second row is the ROI to be resliced
            P = char(mask_path, roi_path);
            spm_reslice(P, struct('mean', false));  % Reslice without estimating a mean image
            % The resliced ROI is saved with prefix 'r'
            [p, n, ext] = fileparts(roi_path);
            roi_resliced_path = fullfile(p, ['r' n ext]);
            if exist(roi_resliced_path, 'file') ~= 2
                warning('Resliced ROI file not found: %s, skipping.', roi_resliced_path);
                continue;
            end
            roi_nii = spm_vol(roi_resliced_path);
            roi_data = spm_read_vols(roi_nii);
            roi_bin = roi_data > 0;
        end
        
        % Calculate voxel overlap and proportion
        overlap = roi_bin & subj_bin;
        overlap_voxels = sum(overlap(:));
        roi_voxels = sum(roi_bin(:));
        if roi_voxels == 0
            prop_overlap = NaN;
        else
            prop_overlap = overlap_voxels / roi_voxels;
        end
        
        % Append results to cell array
        results = [results; {subj, roi_files{r}, overlap_voxels, roi_voxels, prop_overlap}];
    end
end


% Create table from results and save to CSV
T = cell2table(results, 'VariableNames', {'SubjectID', 'MaskVoxelCount'});
disp(T);
csv_file = fullfile(dataroot, 'mask_voxel_count_all_subjects.csv');
writetable(T, csv_file);
fprintf('Results saved to %s\n', csv_file);