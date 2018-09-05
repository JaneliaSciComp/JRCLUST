%--------------------------------------------------------------------------
function S_clu = S_clu_update_(S_clu, clustersToUpdate, P)
    % update cluster waveform and self correlation score

    S0 = get(0, 'UserData');

    % find cluster center
    for iCluster = 1:numel(clustersToUpdate)
        cluster = clustersToUpdate(iCluster);
        spikes = find(S_clu.spikeClusters == cluster);

        S_clu.spikesByCluster{cluster} = spikes;
        S_clu.clusterSites(cluster) = mode(S0.spikeSites(spikes));
        S_clu.nSpikesPerCluster(cluster) = numel(spikes);
    end

    % update mean waveform
    S_clu = clusterMeanWaveforms(S_clu, clustersToUpdate);
    vrSelfCorr_clu = get_diag_(S_clu.mrWavCor);
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, clustersToUpdate);
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, vrSelfCorr_clu);

    for iCluster = 1:numel(clustersToUpdate)
        cluster = clustersToUpdate(iCluster);
        S_clu.mrWavCor(cluster, cluster) = S_clu_self_corr_(S_clu, cluster, S0);
    end

    S_clu = S_clu_position_(S_clu, clustersToUpdate);
    S_clu = S_clu_quality_(S_clu, P, clustersToUpdate);
end % function
