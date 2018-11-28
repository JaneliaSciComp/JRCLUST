%--------------------------------------------------------------------------
function hClust = post_merge_wav_(hClust, spikeData, hCfg)
    % create covariance matrix (mrDist_wav)
    hClust.computeMeanWaveforms(hCfg);
    hClust.mrWavCor = S_clu_wavcor_(hClust, hCfg);

    for iRepeat = 1:hCfg.nPassesMerge %single-pass vs dual-pass correction
        [hClust, nMerges_clu] = S_clu_wavcor_merge_(hClust, hCfg);
        if nMerges_clu < 1, break; end
    end %for
end %func
