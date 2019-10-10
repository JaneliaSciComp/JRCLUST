function phyData = loadPhy(loadPath)
%LOADPHY Load and return Phy-formatted .npy files
phyData = struct();

if exist('readNPY', 'file') ~= 2
    warning('Please make sure you have npy-matlab installed (https://github.com/kwikteam/npy-matlab)');
    return;
end

loadPath_ = jrclust.utils.absPath(loadPath);
if isempty(loadPath_)
    error('Could not find path ''%s''', loadPath);
elseif exist(loadPath, 'dir') ~= 7
    error('''%s'' is not a directory', loadPath);
end

loadPath = loadPath_;
ls = dir(loadPath);
for iFile = 1:numel(ls)
    file = ls(iFile);
    if file.isdir || isempty(regexp(file.name, 'py$', 'once'))
        continue
    end

    [~, fn, ext] = fileparts(file.name);
    switch ext
        case '.py'
            if strcmp(fn, 'params')
                phyData.params = parseParams(fullfile(loadPath, file.name));
            end
            
        case '.npy'
            fn = matlab.lang.makeValidName(fn);
            phyData.(fn) = readNPY(fullfile(loadPath, file.name));
    end
end

phyData.loadPath = loadPath;
end  % function

%% LOCAL FUNCTIONS
function params = parseParams(paramFile)
%PARSEPARAMS Parse a Phy params.py file, returning a struct
params = [];
if exist(paramFile, 'file') ~= 2
    return;
end

params = struct();
keysVals = cellfun(@(line) strsplit(line, '='), jrclust.utils.readLines(paramFile), 'UniformOutput', 0);
for i = 1:numel(keysVals)
    kv = cellfun(@strip, keysVals{i}, 'UniformOutput', 0);
    key = kv{1}; val = kv{2};

    % strip quote marks from strings
    val = regexprep(val, '^r?[''"]', '');  % remove r", r', ", or ' from beginning of string
    val = regexprep(val, '[''"]$', '');  % remove " or ' from end of string

    switch key
        case 'dat_path' % get full path to dat_path
            val_ = jrclust.utils.absPath(val, fileparts(paramFile));
            if isempty(val_)
                warning('File ''%s'' was not found; consider updating params.py with the new path', val);
            end
            val = val_;

        case {'n_channels_dat', 'offset', 'sample_rate', 'n_features_per_channel'}
            val = str2double(val);

        case 'hp_filtered'
            val = strcmp(val, 'True');

        case 'dtype'
            val = strip(val, '''');
    end

    params.(key) = val;
end
end  %function
