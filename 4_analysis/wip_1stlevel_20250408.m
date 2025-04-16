% filepath: /Users/z/Documents/2_STAR_MRI/scripts/mri_preprocessing/4_analysis/wip_1stlevel_20250408
%-----------------------------------------------------------------------
% Design specific parameters
%-----------------------------------------------------------------------

clear; clc
dataroot = '/data/project/STAR/STAR_preprocessed/L_fmriprep/data/derivatives'; % Preprocessed data directory

subs = {
    % 'sub-ID01', 'sub-ID02', 'sub-ID03',... % example subjects
    'sub-ID09', 'sub-ID10'
    % add more subjects as needed
};
nsubs = length(subs);
nsruns = 3;  % Number of runs per session 

cnames{1} = {'OwnCue', 'OwnElaborate', 'OwnRating', 'OtherCue', 'OtherElaborate', 'OtherRating', 'GroundCueOnset'};
ncond = length(cnames{1});
TR = 2; 

hpf = 74; % High-pass filter (sec)
incmoves = 1; % Include movement parameters
modeldur = 1; % Model duration
movefilt = '^sub-.*_ses-mri01_task-TRMrun%d_desc-confounds_timeseries\\.tsv$'; % Movement file filter

mnames_static = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', ...
          'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', ...
          'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1', ...
          'cosine00', 'cosine01', 'cosine02'};
      
failed = {}; % Initialize failed list

for sub = 1:nsubs
    try
        % Define subject's func directory
        subdata = fullfile(dataroot, subs{sub}, 'ses-mri01', 'func');
        cd(subdata);
        if ~exist('analysisext20250408','dir')
            mkdir('analysisext20250408')
        end
        cd('analysisext20250408')
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(subdata, 'analysisext20250408'));
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        for run = 1:nsruns
                        % Check if movement files are there for the current run
            if incmoves == 1
                movefilt_run = sprintf('^sub-.*_ses-mri01_task-TRMrun%d_desc-confounds_timeseries\\.tsv$', run);
                mfiles = spm_select('FPList', subdata, movefilt_run);
                if isempty(mfiles)
                    disp(['No movement files found for run ' num2str(i) ' for subject ' subs{sub}])
                    continue;
                else
                    disp('Movement files found:')
                    disp(mfiles)
                end
                % Filter only the .tsv file in case both JSON and TSV are returned
                mfilesCell = cellstr(mfiles);
                tsv_idx = ~cellfun('isempty', regexp(mfilesCell, '\.tsv$', 'once'));
                if ~any(tsv_idx)
                    error('No TSV movement file found for run %d', i);
                end
                mfiles = char(mfilesCell(tsv_idx));  % use only the TSV file(s)
                
                % Read the movement file (using the first TSV file if multiple)
                movevars = tdfread(mfiles(1,:)); 
                rp = [];
                for p = 1:length(mnames_static)
                    eval(['rp = [rp movevars.' mnames_static{p} '];']);
                end
                % Include spike columns if they exist
                for x = 1:155
                    if isfield(movevars, ['spike_' num2str(x)])
                        eval(['rp = [rp movevars.spike_' num2str(x) '];']);
                    end
                end
                % (Optional) Add compcor columns if present…
                R = rp;
                save(fullfile(subdata, ['movevals_run' num2str(i) '.mat']), 'rp', 'R');
            end
            
            % Update onsets filter (remove extra escape)
            onsetsfilt = sprintf('.*_run%d_pmbh.mat$', i);
            onsetsfile = spm_select('FPList', pwd, onsetsfilt);
            if isempty(onsetsfile)
                error('Onsets/durations file is empty for run %d', i);
            end
            load(onsetsfile)
            
            for k = 1:length(onsets)
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).name = cnames{1}{k};
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).onset = onsets{k};
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).duration = durations{k};
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).tmod = 0;
                if k == 1
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).pmod.name = 'own ratings';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).pmod.param = pmod.OwnRating;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).pmod.poly = 1;
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(k).orth = 1;
            end

            matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).regress = struct('name', {}, 'val', {});
            if incmoves == 1
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi_reg = cellstr(fullfile(subdata, ['movevals_run' num2str(run) '.mat']));
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi_reg = {''};
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).hpf = hpf;
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File',...
            substruct('.','val',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
            substruct('.','val',{2}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'TraumaCue > CtrlCue';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 -1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'TraumaElaborate > CtrlElaborate';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 -1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'TraumaCueElaborate > CtrlCueElaborate';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 1 0 -1 -1];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'Ratings > Cues';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [-1 0 1 -1 0 1];
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
        matlabbatch{5}.spm.stats.con.consess{5}.tcon.name = 'Cue > Elaborate';
        matlabbatch{5}.spm.stats.con.consess{5}.tcon.weights = [1 -1 1 -1];
        matlabbatch{5}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
        matlabbatch{6}.spm.stats.con.consess{6}.tcon.name = 'TraumaCue > Rest';
        matlabbatch{6}.spm.stats.con.consess{6}.tcon.weights = [1 0 0 0 0 0 0];
        matlabbatch{6}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
        matlabbatch{7}.spm.stats.con.consess{7}.tcon.name = 'TraumaElaborate > Rest';
        matlabbatch{7}.spm.stats.con.consess{7}.tcon.weights = [0 1];
        matlabbatch{7}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';
        matlabbatch{7}.spm.stats.con.delete = 0;
        
        cd(fullfile(subdata, 'analysisext20250408'));
        save('firstLevel_batch.mat', 'matlabbatch');
        spm_jobman('run', matlabbatch)
        clear matlabbatch
        
    catch ME
        disp(['Failed processing subject: ' subs{sub}])
        disp(['Error message: ' ME.message])
        failed{end+1} = subs{sub};
        continue;
    end
end

if ~isempty(failed)
    disp('The following subjects failed processing:')
    disp(failed)
else
    disp('All subjects processed successfully!')
end