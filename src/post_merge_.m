%--------------------------------------------------------------------------
function [S_clu, S0] = post_merge_(S_clu, P, fPostCluster)
    % waveform based merging. find clusters within maxSite
    % also removes duplicate spikes
    if nargin<3, fPostCluster=1; end

    nRepeat_merge = get_set_(P, 'nRepeat_merge', 10);
    % fMerge_pv = 0;
    fClean_clu = 1;

    if fPostCluster, S_clu = postCluster_(S_clu, P); end

    % Add S_clu fields
    % S0 = get(0, 'UserData');
    S_clu = rmfield_(S_clu, 'tmrWav_clu', 'vrPosX_clu', 'vrPosY_clu', 'vrVpp_clu', ...
    'vrVpp_uv_clu', 'vrVmin_uv_clu', 'vrSnr_clu', 'vnSite_clu', ...
    'vrIsoDist_clu', 'vrLRatio_clu', 'vrIsiRatio_clu');
    S_clu = S_clu_refresh_(S_clu);
    S_clu = S_clu_sort_(S_clu, 'viSite_clu');
    S_clu = rmfield_(S_clu, 'csNote_clu');

    if fClean_clu, S_clu = S_clu_cleanup_(S_clu, P); end
    S_clu = post_merge_wav_(S_clu, nRepeat_merge, P);
    S_clu = S_clu_refresh_(S_clu);
    S_clu = S_clu_sort_(S_clu, 'viSite_clu');
    S_clu = S_clu_update_wav_(S_clu, P);

    % set diagonal element
    [S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, [], S0));
    S_clu.P = P;
    S_clu = S_clu_position_(S_clu);
    S_clu.csNote_clu = cell(S_clu.nClu, 1); %reset note
    S_clu = S_clu_quality_(S_clu, P);
    [S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
end %func
