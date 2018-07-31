%--------------------------------------------------------------------------
function S_ksort = kilosort(vcFile_prm)
    % add actual KiloSort repo to path
    [dirname, ~] = fileparts(fullfile(mfilename('fullpath')));
    if exist(fullfile(dirname, 'KiloSort'), 'dir') == 7
        addpath(genpath(fullfile(dirname, 'KiloSort')));
    elseif ~exist('preprocessData', 'file') ...
        || ~exist('fitTemplates', 'file') ...
        || ~exist('fullMPMU', 'file')
        error('KiloSort functions not found');
    end

    % Run Kilosort
    if ischar(vcFile_prm)
        P = loadParams(vcFile_prm);
        if isempty(P), return; end
    else
        P = vcFile_prm;
        vcFile_prm = P.prmFile;
    end

    fprintf('Running KiloSort on %s\n', vcFile_prm);
    
    fSave_phy = get_set_(P, 'fSave_phy', 0);
    fMerge_post = get_set_(P, 'fMerge_post', 0);
    
    runtime_ksort = tic;

    % Run KiloSort
    [fpath, ~, ~] = fileparts(vcFile_prm);
    if isempty(fpath)
        fpath = pwd;
    end

    ops = kilosort_config(P); % get config
    S_chanMap = kilosort_chanMap(P); %make channel map

    [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
    rez                = fitTemplates(rez, DATA, uproj); % fit templates iteratively
    rez                = fullMPMU(rez, DATA); % extract final spike times (overlapping extraction)

    if fMerge_post && exist('merge_posthoc2', 'file')
        rez = merge_posthoc2(rez);
    elseif fMerge_post
        warning('merge_posthoc2 not found; skipping');
    end

    % save python results file for Phy
    if fSave_phy
        if exist(fullfile(dirname, 'npy-matlab'), 'dir')
            addpath(fullfile(dirname, 'npy-matlab'));
        end

        if ~exist('rezToPhy', 'file')
            warning('rezToPhy not found; skipping');
        else

            try
                rezToPhy(rez, fpath); %path to npy2mat needed
            catch
                disperr_();
            end
        end
    end

    runtime_ksort = toc(runtime_ksort);
    fprintf('\tKiloSort took %0.1fs for %s\n', runtime_ksort, P.prmFile);

    % output KiloSort result
    S_ksort = struct('rez', rez, 'P', P, 'runtime_ksort', runtime_ksort);
    struct_save_(S_ksort, strrep(vcFile_prm, '.prm', '_ksort.mat'), 1);
end %func
