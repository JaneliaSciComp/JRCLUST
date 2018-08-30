%--------------------------------------------------------------------------
function S_clu = mergePair(S_clu, P, cluster1, cluster2)
    % merge a pair of clusters

    % ensure cluster1 < cluster2
    if cluster1 > cluster2
        tmp = cluster2;
        cluster2 = cluster1;
        cluster1 = tmp;
    end

    % cluster1 takes all spikes; cluster2 gets nothing
    S_clu.nSpikesPerCluster(cluster1) = sum(S_clu.nSpikesPerCluster([cluster1 cluster2]));
    S_clu.nSpikesPerCluster(cluster2) = 0;

    % update spike assignments
    S_clu.spikeClusters(S_clu.spikeClusters == cluster2) = cluster1;
    S_clu.spikesByCluster{cluster1} = find(S_clu.spikeClusters == cluster1);
    S_clu.spikesByCluster{cluster2} = [];

    try % TODO: this shouldn't need a try statement
        S_clu.clusterNotes{cluster1} = '';
        S_clu.clusterNotes{cluster2} = '';
    catch
    end

    S_clu = removeSpikesInRefracPeriod(S_clu, P, cluster1);
    S_clu = S_clu_update_(S_clu, cluster1, P);
    S_clu = delete_clu_(S_clu, cluster2);
    S_clu = S_clu_remove_empty_(S_clu);
end
