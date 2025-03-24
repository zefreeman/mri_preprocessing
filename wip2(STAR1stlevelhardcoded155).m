%-----------------------------------------------------------------------
% Design specific parameters
%-----------------------------------------------------------------------

clear

subs = {'sub-ID61'}; % List of subjects
nsubs = length(subs);
nruns = 3; % Number of runs per session

sessname = {'TraumaReminder'};
eventroot = '/data/project/STAR/STAR_preprocessed/N_fmriprep/data/derivatives'; % Onset and durs directory
dataroot = '/data/project/STAR/STAR_preprocessed/N_fmriprep/data/derivatives'; % Preprocessed data directory

cnames{1} = {'OwnCue', 'OwnElaborate', 'OwnRating', 'OtherCue', 'OtherElaborate', 'OtherRating', 'GroundCueOnset'};
ncond = length(cnames{1});
TR = 2; % RT (in seconds)

hpf = 74; % High-pass filter (in seconds)
incmoves = 1; % Include movement parameters
modeldur = 1; % Model duration
imgfilt = '^sub-.*_ses-mri01_task-TRMrun%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii$'; % Image filter
movefilt = '^sub-.*_ses-mri01_task-TRMrun%d_desc-confounds_timeseries.tsv$'; % Movement file filter

% Columns to include as regressors
mnames = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', ...
          'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', ...
          'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1', ...
          'a_comp_cor_00', 'a_comp_cor_01', 'a_comp_cor_02', 'a_comp_cor_03', ...
          'a_comp_cor_04', 'a_comp_cor_05', 'a_comp_cor_06', 'a_comp_cor_07', ...
          'white_matter', 'csf', 'cosine00', 'cosine01', 'cosine02'};

%-----------------------------------------------------------------------
% Design setup
%-----------------------------------------------------------------------

% Basis functions and timing parameters
%---------------------------------------------------------------------------
xBF.name     = 'hrf';
xBF.length   = 32.2;              % Length in seconds
xBF.order    = 1;                   % Order of basis set
xBF.T        = 40;                    % Microtime resolution (number of time bins per scan)
xBF.T0       = 20;                   % Microtime onset (time bin at which responses are modeled)
xBF.UNITS    = 'secs';           % Units of onsets and durations
xBF.Volterra = 1;                  % Order of convolution (1 = no derivatives)
%!!!% think there is an issue here with X - it should eventually have the
%number of rows of runs 1,2,3 - 
failed = {}; % Initialize failed list

for sub = 1:nsubs
    try
        clear SPM
        disp(['Processing subject: ' subs{sub}])
        SPM.xY.RT = TR;
        SPM.xGX.iGXcalc = 'None';
        SPM.xVi.form = 'AR(1)';
        SPM.xBF = xBF;
        csub = subs{sub};

        subdata = fullfile(dataroot, csub, 'ses-mri01'); % Subject data directory
        anadir = fullfile(subdata, 'func'); % Analysis directory
    % if exist(anadir, 'dir') ~= 7; mkdir(subdata, 'func'); end % Create
    % func directory if it doesn't exist - although why wouldn't it

        cd(anadir);
        tc = 0; % Initialise time course counter - breaks if not here, unsure on why given there are other loops
        allfiles = ''; % Initialise list of al functional files

        for run = 1:nruns
            evorder = cnames{1}; % Assuming one session
            clear ffiles;

            sessdata = fullfile(subdata, 'func');

            % Select functional image files for the current run
            imgfilt_run = sprintf(imgfilt, run);
            files = spm_select('List', sessdata, imgfilt_run);
            if isempty(files)
                disp(['No functional image files found for run ' num2str(run)])
                continue;
            else
                disp('Functional image files found:')
                disp(files)
            end
            
            % Prepare movement files for the current run
            if incmoves == 1
                movefilt_run = sprintf(movefilt, run);
                mfiles = spm_select('List', sessdata, movefilt_run);
                if isempty(mfiles)
                    disp(['No movement files found for run ' num2str(run)])
                    continue;
                else
                    disp('Movement files found:')
                    disp(mfiles)
                end
                
                % Convert multiple file names into cell array
                mfiles = cellstr(mfiles);
                movefiles = fullfile(sessdata, mfiles);
                
                % Initialize moves array
                moves = [];
                
                % Read each movement file and extract required columns
                for i = 1:length(movefiles)
                    movefile = movefiles{i};
                    movevars = readtable(movefile, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', true);

                    for col = 1:width(movevars)
                           % Replace 'n/a' with 0 and convert column by
                           % column to numeric
                           if ~isnumeric(movevars{:, col});
                               cellArray =movevars{:, col};
                               cellArray = strrep(cellArray, 'n/a', '0');
                               % Convert to numeric here
                               cellArray = str2double(cellArray);
                               movevars.(col) = cellArray;
                           end
                    end
                    
                    % Check and extract required columns
                    movevars = movevars(:, mnames);
                    if any(varfun(@iscell, movevars,'OutputFormat','uniform'))
                        error(['Inconsistent column types in movement file: ' movefile])
                    end
                    moves = [moves table2array(movevars)];
                end
            end

            for f = 1:size(files, 1)
                ffiles{f} = fullfile(sessdata, files(f, :));
            end
            allfiles = [allfiles; ffiles']; % Concatenate file paths

            % Load the .mat file containing onsets and durations for the current run
            onsdur_pattern = sprintf('.*onsets_durations_run%d\\.mat$', run);
            onsdurfile = spm_select('FPList', sessdata, onsdur_pattern);
            if isempty(onsdurfile)
                error(['onsets_durations file not found for run ' num2str(run) ' in session data directory: ' sessdata])
            else
                disp(['Onsets and durations file found for run ' num2str(run) ':'])
                disp(onsdurfile)
            end
            load(onsdurfile, 'onsets', 'durations', 'names'); 
            
            tc = tc + 1; % Increment time course counter - NOT SURE WHAT WENT WRONG BUT THIS BROKE FOR A WHILE - IS SOURCE OF OVERWRITE?
            for c = 1:ncond
                if length(onsets{c}) ~= length(durations{c})
                    error(['Inconsistent onset and duration lengths for condition: ' cnames{1}{c}])
                end
                SPM.Sess(tc).U(c) = struct('ons', onsets{c}, 'dur', durations{c}, 'name', {names(c)}, 'P', struct('name', 'none'));
            end

            SPM.nscan(tc) = 155; % Number of scans - volumes - trying '155' instead of = size(files,1) - THIS IS NOT A LONG TERM FIX 
            SPM.xX.K(tc).HParam = hpf; % High pass filter parameter

            if incmoves == 1
                SPM.Sess(tc).C.C = moves;
                SPM.Sess(tc).C.name = mnames;
            else
                SPM.Sess(tc).C.C = [];
                SPM.Sess(tc).C.name = {};
            end
        end
        SPM.xY.P = char(allfiles); % List of all functional files

        cd(anadir)
        SPMdes = spm_fmri_spm_ui(SPM); % Design specification
        spm_unlink(fullfile('.', 'maskFirstLevel.img')); % Avoid overwrite dialog - THIS DOESN'T WORK 
        SPMest = spm_spm(SPMdes); % Model estimation
    catch ME
        disp(['Failed processing subject: ' subs{sub}])
        disp(['Error message: ' ME.message])
        failed{end+1} = subs{sub}; % Record failed subject
        % Optionally, log stack trace for more detailed debugging
        disp(getReport(ME, 'extended'));
    end
end

disp('Failed subjects:');
disp(failed);

