%--------------------------------------------------------------------------
function S_ksort = kilosort(vcFile_prm)
    % add actual KiloSort repo to path
    [dirname, ~] = fileparts(fullfile(mfilename('fullpath')));
    addpath(genpath(fullfile(dirname, 'KiloSort')));

    % Run Kilosort
    if ischar(vcFile_prm)
        fprintf('Running kilosort on %s\n', vcFile_prm);
        P = loadParam_(vcFile_prm);
        if isempty(P), return; end
    else
        P = vcFile_prm;
        vcFile_prm = P.vcFile_prm;
        fprintf('Running kilosort on %s\n', vcFile_prm);
    end
    
    fSave_phy = get_set_(P, 'fSave_phy', 0);
    fMerge_post = get_set_(P, 'fMerge_post', 0);
    
    runtime_ksort = tic;

    % Run kilosort
    [fpath, ~, ~] = fileparts(vcFile_prm);
    ops = kilosort_config(P); % get config
    S_chanMap = kilosort_chanMap(P); %make channel map

    [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
    rez                = fitTemplates(rez, DATA, uproj); % fit templates iteratively
    rez                = fullMPMU(rez, DATA); % extract final spike times (overlapping extraction)

    if fMerge_post
        rez = merge_posthoc2(rez);
    end

    % save python results file for Phy
    if fSave_phy
        addpath(fullfile(dirname, 'npy-matlab'));
        try
            rezToPhy(rez, fpath); %path to npy2mat needed
        catch
            disperr_();
        end
    end

    runtime_ksort = toc(runtime_ksort);
    fprintf('\tkilosort took %0.1fs for %s\n', runtime_ksort, P.vcFile_prm);

    % output kilosort result
    S_ksort = struct('rez', rez, 'P', P, 'runtime_ksort', runtime_ksort);
    struct_save_(S_ksort, strrep(vcFile_prm, '.prm', '_ksort.mat'), 1);
end %func
