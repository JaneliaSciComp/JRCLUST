%--------------------------------------------------------------------------
function S0 = sort_(P, fLoad)
    % Extract feature and sort
    runtime_sort = tic;
    if nargin<2, fLoad = 1; end
    if fLoad
        [S0, P] = load_cached_(P);
    else
        S0 = get(0, 'UserData');
    end

    % Sort and save
    S_clu = fet2clu_(S0, P);
    if get_set_(P, 'fCorrect_overlap', 0) % correct waveforms and features after correcting clusters
        S_clu = sort_overlap_(S0, S_clu, P);
    end
    [S_clu, S0] = S_clu_commit_(S_clu, 'sort_');
    % S0 = set0_(P); %, dimm_fet, cvrTime_site, cvrVpp_site, cmrFet_site, P);

    % measure time
    runtime_sort = toc(runtime_sort);
    fprintf('Sorting took %0.1fs for %s\n', runtime_sort, P.vcFile_prm);
    S0 = set0_(runtime_sort, P);
    S0 = clear_log_(S0);

    save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
end %func
