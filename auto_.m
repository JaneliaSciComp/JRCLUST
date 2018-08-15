%--------------------------------------------------------------------------
function auto_(P)
    if getOr(P, 'fRepeat_clu', 0)
        sort_(P);
        describe_(P.paramFile);
        return;
    end
    [S0, P] = load_cached_(P); % load cached data or from file if exists
    S_clu = get_(S0, 'S_clu');
    S_clu.P = P;

    if isempty(S_clu)
        fprintf(2, 'You must sort first by running "jrc sort".\n');
        return;
    end
    [S_clu, S0] = post_merge_(S_clu, P);
    S0 = clear_log_(S0);
    save0_(strrep(P.paramFile, '.prm', '_jrc.mat'));
end % func
