%--------------------------------------------------------------------------
function flag = is_sorted_(P)
    % return true if already detected. .spkwav file must exist
    S0 = load0_(strrep(P.prmFile, '.prm', '_jrc.mat'));
    S_clu = get_(S0, 'S_clu');
    flag = ~isempty(S_clu);
end %func
