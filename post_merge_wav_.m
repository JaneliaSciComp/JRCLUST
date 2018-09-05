%--------------------------------------------------------------------------
function S_clu = post_merge_wav_(S_clu, nRepeat_merge, P)
    % S0 = get(0, 'UserData');

    % create covariance matrix (mrDist_wav)
    S_clu = clusterMeanWaveforms(S_clu);
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P);

    for iRepeat = 1:nRepeat_merge %single-pass vs dual-pass correction
        [S_clu, nMerges_clu] = S_clu_wavcor_merge_(S_clu, P);
        if nMerges_clu < 1, break; end
    end %for
end % function
