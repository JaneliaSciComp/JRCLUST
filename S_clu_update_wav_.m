%--------------------------------------------------------------------------
function S_clu = S_clu_update_wav_(S_clu, P)
    if nargin<2, P = get0_('P'); end
    % if nargin<2, S0 = get(0, 'UserData'); end

    S_clu = S_clu_wav_(S_clu);
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P);
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu));
end %func
