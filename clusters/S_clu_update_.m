%--------------------------------------------------------------------------
function S_clu = S_clu_update_(S_clu, clustersToUpdate, P)
    % update cluster waveform and self correlation score
    % mrWav not needed
    S0 = get(0, 'UserData');

    % find cluster center
    for iClu = 1:numel(clustersToUpdate)
        cluster = clustersToUpdate(iClu);
        spikesThisCluster = find(S_clu.spikeClusters == cluster);
        % why do all this twice? TODO: investigate
        S_clu.spikesByCluster{cluster} = spikesThisCluster;
        S_clu.clusterSites(cluster) = mode(S0.spikeSites(spikesThisCluster));
        S_clu.nSpikesPerCluster(cluster) = numel(spikesThisCluster);
    end

    % update mean waveform
    S_clu = S_clu_wav_(S_clu, clustersToUpdate);
    
    vrSelfCorr_clu = get_diag_(S_clu.mrWavCor);
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, clustersToUpdate);
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, vrSelfCorr_clu);
    for iClu = 1:numel(clustersToUpdate)
        cluster = clustersToUpdate(iClu);
        S_clu.mrWavCor(cluster,cluster) = S_clu_self_corr_(S_clu, cluster, S0);
    end
    S_clu = S_clu_position_(S_clu, clustersToUpdate);
    S_clu = S_clu_quality_(S_clu, P, clustersToUpdate);
    % [S_clu, S0] = S_clu_commit_(S_clu, 'S_clu_update_');
end % function
