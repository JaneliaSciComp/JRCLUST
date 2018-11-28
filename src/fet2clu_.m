%--------------------------------------------------------------------------
function S_clu = fet2clu_(S0, P)
    % can process different shanks separately
    fprintf('Clustering\n');
    S_clu = cluster_spacetime_(S0, P);
    S_clu = postCluster_(S_clu, P);

    if get_set_(P, 'fRepeat_clu', 0), S_clu = S_clu_reclust_(S_clu, S0, P); end

    spikeData = struct('spikeTimes', S0.viTime_spk, ...
                    'spikeSites', S0.viSite_spk, ...
                    'spikeSites2', S0.viSite2_spk, ...
                    'spikePositions', S0.mrPos_spk);
    S_clu = jrclust.clustering.autoMerge(S_clu, spikeData, P, 0);

    S_clu.viClu_auto = S_clu.viClu;
    fprintf('\tClustering took %0.1f s\n', S_clu.t_runtime);
end %func
