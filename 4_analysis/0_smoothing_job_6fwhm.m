%-----------------------------------------------------------------------
% Job saved on 13-Aug-2024 14:43:02 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

clear; clc

dataroot = '/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/derivatives'; % Preprocessed data directory

subs = {'sub-ID01', 'sub-ID02', 'sub-ID03',... %'sub-ID04', ... % London pts 
  'sub-ID05', 'sub-ID06', 'sub-ID07', 'sub-ID08',...
  'sub-ID09', 'sub-ID10',...%, 'sub-ID11', 'sub-ID12'
  'sub-ID13', 'sub-ID14', 'sub-ID15', 'sub-ID16', ...
  'sub-ID17',...% 'sub-ID18', 
  'sub-ID19',... %'sub-ID20', ...
  'sub-ID21',... %'sub-ID22', 
  'sub-ID23',... %'sub-ID24', ...
  'sub-ID25', 'sub-ID26', 'sub-ID62',... %'sub-ID63', ...
  'sub-ID64', 'sub-ID65', 'sub-ID66', 'sub-ID67', ...
  'sub-ID68', 'sub-ID69', 'sub-ID70', 'sub-ID71', ...
  'sub-ID72',... %'sub-ID73', 
  'sub-ID74'    
% 'sub-ID47', 'sub-ID48', 'sub-ID49', 'sub-ID50', ... % Newcastle pts
% 'sub-ID51', 'sub-ID52', 'sub-ID53', 'sub-ID54', ...
% 'sub-ID55', 'sub-ID56', 'sub-ID57', 'sub-ID58', ...
% 'sub-ID59', 'sub-ID60', 'sub-ID61', 'sub-ID62' 
%             'sub-ID27', 'sub-ID31', 'sub-ID32',... Manchester pts
%             'sub-ID34', 'sub-ID35', 'sub-ID39',...
%             'sub-ID41', 'sub-ID44', 'sub-ID45',...
%             'sub-ID37', 'sub-ID42', 'sub-ID46'
            };   % List of subjects
nsubs = length(subs);
nsruns = 3; % Number of runs per session

failed = {}; % Initialize failedlist

for sub = 1:nsubs
    try
        subdata = [dataroot filesep subs{sub} '/ses-mri02/func/'];

        for i = 1:nsruns

         % Define the filter to pick the preprocessed task images for smoothing
         imgfilt = ['^sub-*.*_ses-mri02_task-TRMrun' num2str(i) '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii$']; % Image filter
         files = spm_select('ExtFPList', subdata, imgfilt);

            % Initialise the batch 
            matlabbatch{1}.spm.spatial.smooth.data = cellstr(files);
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's6';

            % Execute the batch
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
            
        end
        
      catch ME
      % Log the subject that caused the error
      disp(['Failed processing subject: ' subs{sub}])  % I don't think this is working properly - should  have been some failed out of M
      disp(['Error message: ' ME.message])
      failed{end+1} = subs{sub}; % Add the failed subject to the list
      continue; % Move on to the next subject
        
    end
    
end