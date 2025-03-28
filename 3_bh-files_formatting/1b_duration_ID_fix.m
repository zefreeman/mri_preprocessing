% Define the folder where the .mat files are located
folder = '/data/project/STAR/STAR_preprocessed/L_OnsDur/';  % Replace with the actual folder path where your .mat files are located

% Define the file names to be renamed and their corresponding new names
fileNames = {
    'P040093', 'sub-ID19';
    'P010095', 'sub-ID04';
    'P020139', 'sub-ID27';
    'P030124', 'sub-ID18';
    'P020149', 'sub-ID28';
    'P040130', 'sub-ID20';
    'P020156', 'sub-ID29';
    'P040155', 'sub-ID21';
    'P020178', 'sub-ID31';
    'P010183', 'sub-ID01';
    'P040181', 'sub-ID03';
    'P020172', 'sub-ID30';
    'P020179', 'sub-ID32';
    'P020189', 'sub-ID33';
    'P010200', 'sub-ID06';
    'P010201', 'sub-ID02';
    'P020190', 'sub-ID34';
    'P010209', 'sub-ID07';
    'P020202', 'sub-ID35';
    'P030208', 'sub-ID47';
    'P030215', 'sub-ID48';
    'P040225', 'sub-ID64';
    'P050196', 'sub-ID71';
    'P040246', 'sub-ID67';
    'P010242', 'sub-ID08';
    'P020228', 'sub-ID36';
    'P010162', 'sub-ID05';
    'P040230', 'sub-ID65';
    'P050237', 'sub-ID72';
    'P020252', 'sub-ID37';
    'P040270', 'sub-ID68';
    'P040231', 'sub-ID66';
    'P030285', 'sub-ID49';
    'P010314', 'sub-ID09';
    'P040310', 'sub-ID69';
    'P040320', 'sub-ID70';
    'P010337', 'sub-ID11';
    'P010327', 'sub-ID10';
    'P030335', 'sub-ID52';
    'P020336', 'sub-ID39';
    'P030287', 'sub-ID50';
    'P030342', 'sub-ID53';
    'P030346', 'sub-ID54';
    'P020331', 'sub-ID38';
    'P050328', 'sub-ID73';
    'P020339', 'sub-ID40';
    'P010345', 'sub-ID12';
    'P050355', 'sub-ID74';
    'P050352', 'sub-ID23';
    'P020358', 'sub-ID41';
    'P020359', 'sub-ID42';
    'P040372', 'sub-ID22';
    'P030312', 'sub-ID51';
    'P050368', 'sub-ID24';
    'P020369', 'sub-ID43';
    'P020394', 'sub-ID45';
    'P050377', 'sub-ID25';
    'P010380', 'sub-ID13';
    'P010386', 'sub-ID14';
    'P020392', 'sub-ID44';
    'P030398', 'sub-ID56';
    'P010399', 'sub-ID15';
    'P010405', 'sub-ID16';
    'P050417', 'sub-ID26';
    'P010420', 'sub-ID17';
    'P030421', 'sub-ID57';
    'P010426', 'sub-ID63';
    'P020431', 'sub-ID46';
    'P030397', 'sub-ID55';
    'P030434', 'sub-ID61';
    'P010425', 'sub-ID62';
    'P030432', 'sub-ID59';
    'P030433', 'sub-ID60'
};

% Loop through each file name pair
for k = 1:size(fileNames, 1)
    oldFileName = fileNames{k, 1};
    newFileName = fileNames{k, 2};

    % Generate the regular expression pattern to match the old file name
    pattern = ['(?<=' oldFileName ').*'];

    % List all .mat files in the folder
    files = dir(fullfile(folder, ['*' oldFileName '*']));

    % Loop through each .mat file
    for i = 1:length(files)
        % Check if the file name matches the pattern
        if ~isempty(regexp(files(i).name, pattern, 'once'))
            % Remove 'STAR' if present and replace oldFileName with newFileName
            newName = strrep(files(i).name, 'STAR', '');
            newName = regexprep(newName, [oldFileName '.*'], [newFileName '_task-TRM_$0']);

            % Append _bh before the extension
            [~, baseFileName, ext] = fileparts(newName);
            newFullFileName = fullfile(folder, [baseFileName '_bh' ext]);

            % Construct the full old file name and new file name
            oldFullFileName = fullfile(folder, files(i).name);

            % Rename the file by moving it to the new file name
            movefile(oldFullFileName, newFullFileName);
        end
    end
end
