%-----------------------------------------------------------------------
% Ze Freeman 202407 - 202504
%-----------------------------------------------------------------------
%copy of Xruns_20241218_1stlevel_job.m which was the version with the
%original number of contrasts
 
clear; clc
dataroot = '/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/derivatives'; % Preprocessed data directory

subs = {
'sub-ID01', 'sub-ID02', 'sub-ID03', 'sub-ID04', ... % London pts 
'sub-ID05', 'sub-ID06', 'sub-ID07', 'sub-ID08', ...
'sub-ID09', 'sub-ID10', 'sub-ID11', 'sub-ID12', ...
'sub-ID13', 'sub-ID14', 'sub-ID15', 'sub-ID16', ...
'sub-ID17', 'sub-ID18', 'sub-ID19', 'sub-ID20', ...
'sub-ID21', 'sub-ID22', 'sub-ID23', 'sub-ID24', ...
'sub-ID25', 'sub-ID26', 'sub-ID62', 'sub-ID63', ...
'sub-ID64', 'sub-ID65', 'sub-ID67', ...
'sub-ID68', 'sub-ID69', 'sub-ID70', 'sub-ID71', ...
'sub-ID72', 'sub-ID73', 'sub-ID74'
% The following subjects failed processing:
%     'sub-ID19'    'sub-ID25'    'sub-ID63'    'sub-ID67'
% 'sub-ID47', 'sub-ID48', 'sub-ID49', 'sub-ID50', ... % Newcastle pts
% 'sub-ID51', 'sub-ID52', 'sub-ID53', 'sub-ID54', ...
% 'sub-ID56', 'sub-ID57', 'sub-ID58', ...
% 'sub-ID59', 'sub-ID60', 'sub-ID61' 
%'sub-ID60'
% The following subjects failed processing:
%     'sub-ID51'    'sub-ID60'
% 'sub-ID27', 'sub-ID28', 'sub-ID29', 'sub-ID30', ... % Manchester pts
% 'sub-ID33', 'sub-ID35', 'sub-ID36', 'sub-ID37', ...
% 'sub-ID38', 'sub-ID39', 'sub-ID40', 'sub-ID42', ...
% 'sub-ID43', 'sub-ID45', 'sub-ID46'
% The following subjects failed processing:
%     'sub-ID30'    'sub-ID33'    'sub-ID35'    'sub-ID36'    'sub-ID38'    'sub-ID39'
};    
nsubs = length(subs);
nsruns = 3; % Number of runs per session - maximum
%run_files = dir(fullfile(subdata, 'sub-*ses-mri01_task-TRMrun*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'));
%nsruns = numel(run_files); % Automatically detects the number of runs

cnames{1} = {'OwnCue', 'OwnElaborate', 'OwnRating', 'OtherCue', 'OtherElaborate', 'OtherRating', 'GroundCueOnset'};
ncond = length(cnames{1});
TR = 2; 

hpf = 74; % High-pass filter (in seconds)
incmoves = 1; % Include movement parameters
modeldur = 1; % Model duration
movefilt = '^sub-*.*_ses-mri02_task-TRMrun%d_desc-confounds_timeseries.tsv$'; % Movement file filter

% Static columns to include as movement variables - 
% dynamic compcor and others are below
mnames_static = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', ...
          'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', ...
          'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1', ...
          'cosine00', 'cosine01', 'cosine02'};
      
failed = {}; % Initialise failedlist

for sub = 1:nsubs
    try
        % Loop over subjects
            subdata = [dataroot filesep subs{sub} '/ses-mri02/func/'];
            cd(subdata);
            mkdir T2analysis_20250415
            cd T2analysis_20250415
            delete SPM.mat 
            cd(subdata)
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr([subdata 'T2analysis_20250415']); %% 
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
                  
            
             runs = 0;
             for i = 1:nsruns
                 
            % Grab - data, onsets and movement
             imgfilt = ['^s6sub-*.*_ses-mri02_task-TRMrun' num2str(i) '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii$']; % s6 for smoothed 6mm
             files = spm_select('ExtFPList', subdata, imgfilt);
             onsetsfilt = sprintf('.*run%d_pmmi_bh\\.mat$', i);  % 
             onsetsfile = spm_select('FPList', pwd, onsetsfilt);
             load(onsetsfile)
                    
                 if  length(unique(onsets{1})) > 1 %placeholder 'onsets' - previously pmod but 
                     runs = runs+1;
                % Check if movement files are there for the current run
                        if incmoves == 1
                            movefilt_run = sprintf(movefilt, i);
                            mfiles = spm_select('FPList', subdata, movefilt_run);
                            if isempty(mfiles)
                                disp(['No movement files found for run ' num2str(run)])
                                continue;
                            else
                                disp('Movement files found:')
                                disp(mfiles)
                            end
                        end 
                            % Convert multiple file names into cell array
                            % Initialize moves array
                            moves = []; % Needed?

                            % Read each movement file and extract static columns
                                movevars = tdfread(mfiles); 
                                rp = [];
                                present_compcor = [];
                                for p  = 1:length(mnames_static)
                                    eval(['rp = [rp movevars.' mnames_static{p} '];']) ;
                                end

                                % Iterates through all possible spike
                                % number columns and adds if present
                                for x = 1:155
                                        if isfield(movevars,['spike_' num2str(x)])
                                            eval(['rp =[rp movevars.spike_' num2str(x) '];' ]);
                                        end
                                end
                                
                                % Adds up to 5 compcor if they are present
                                for name = {'c_comp_cor_00', 'c_comp_cor_01', ...
                                                                    'c_comp_cor_02', 'c_comp_cor_03', ...
                                                                    'c_comp_cor_04', 'w_comp_cor_00', ... 
                                                                    'w_comp_cor_01', 'w_comp_cor_02', ...
                                                                    'w_comp_cor_03', 'c_comp_cor_04'}
                                         name = string(name);
                                         if isfield(movevars, name)
                                             eval(['rp =[rp movevars.(name)]']);
                                             present_compcor = [present_compcor name];
                                         end
                                end

                                % Make the movement values file per run
                                R=rp; % Fixed a problem 20240813 
                                eval(['save movevals_run' num2str(i) '.mat rp R']); 

            matlabbatch{1}.spm.stats.fmri_spec.sess(runs).scans = cellstr(files);
            
            for k = 1:length(onsets)
                matlabbatch{1}.spm.stats.fmri_spec.sess(runs).cond(k).name =names{k}; %% here? "
                matlabbatch{1}.spm.stats.fmri_spec.sess(runs).cond(k).onset =onsets{k}; %% here? 
                matlabbatch{1}.spm.stats.fmri_spec.sess(runs).cond(k).duration =durations{k};
                matlabbatch{1}.spm.stats.fmri_spec.sess(runs).cond(k).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(runs).cond(k).orth = 1;
                
            end
        matlabbatch{1}.spm.stats.fmri_spec.sess(runs).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(runs).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(runs).multi_reg = cellstr([pwd filesep 'movevals_run' num2str(i) '.mat']) ;
        matlabbatch{1}.spm.stats.fmri_spec.sess(runs).hpf = 74;
        end
                 end
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'OwnCue > OtherCue';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 -1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc'; 
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'OwnElaborate > OtherElaborate';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 -1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'OwnCueElaborate > OtherCueElaborate';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 1 0 -1 -1];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'OwnCueElaborate > Remaining';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 1];
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc'; 
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'OtherCueElaborate > Remaining';
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 -1 -1];
        matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
        % additional checks below
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Ratings > Cues';
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [-1 0 1 -1 0 1];
        matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Cue > Elaborate';
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [1 -1 1 -1];
        matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'TraumaCue > Rest';
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1 0 0 0 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'TraumaElaborate > Rest';
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 1];
        matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.delete = 0;
        
        cd T2analysis_20250415 
        save firstLevel_batch matlabbatch
        spm_jobman('run', matlabbatch)
        clear matlabbatch
        
        catch ME
        % Log the subject that caused the error
        disp(['Failed processing subject: ' subs{sub}])
        disp(['Error message: ' ME.message])
        failed{end+1} = subs{sub}; % Add the failed subject to the list
        continue; % Move on to the next subject
    end
end

if ~isempty(failed)
    disp('The following subjects failed processing:')
    disp(failed)
else
    disp('All subjects processed successfully!')
end

