%--------------------------------------------------------------------------
function S_ksort = kilosort_(vcFile_prm)
    % Run Kilosort
    fSavePhy = 1;
    fMerge_post = 1;
    if ischar(vcFile_prm)
        fprintf('Running kilosort on %s\n', vcFile_prm);
        P = loadParam_(vcFile_prm);
        if isempty(P), return; end
    else
        P = vcFile_prm;
        vcFile_prm = P.vcFile_prm;
        fprintf('Running kilosort on %s\n', vcFile_prm);
    end
    runtime_ksort = tic;

    % Run kilosort
    [fpath, ~, ~] = fileparts(vcFile_prm);
    ops = kilosort('config', P); %get config
    S_chanMap = kilosort('chanMap', P); %make channel map

    [rez, DATA, uproj] = kilosort('preprocessData', ops); % preprocess data and extract spikes for initialization
    rez                = kilosort('fitTemplates', rez, DATA, uproj); % fit templates iteratively
    rez                = kilosort('fullMPMU', rez, DATA); % extract final spike times (overlapping extraction)

    if fMerge_post
        rez = kilosort('merge_posthoc2', rez); %ask whether to merge or not
    end

    % save python results file for Phy
    if fSavePhy
        try
            kilosort('rezToPhy', rez, fpath); %path to npy2mat needed
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
