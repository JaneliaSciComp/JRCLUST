%--------------------------------------------------------------------------
function S_clu = merge_clu_pair_(S_clu, iClu1, iClu2)
    % if iClu1>iClu2, [iClu1, iClu2] = swap(iClu1, iClu2); end

    % update nSpikesPerCluster, spikeClusters, clusterSites. move iClu2 to iClu1
    n1 = S_clu.nSpikesPerCluster(iClu1);
    n2 = S_clu.nSpikesPerCluster(iClu2);
    S_clu.nSpikesPerCluster(iClu1) = n1 + n2;
    S_clu.nSpikesPerCluster(iClu2) = 0;
    S_clu.spikeClusters(S_clu.spikeClusters == iClu2) = iClu1;
    S_clu.spikesByCluster{iClu1} = find(S_clu.spikeClusters == iClu1);
    S_clu.spikesByCluster{iClu2} = [];
    try
        S_clu.clusterNotes{iClu1} = '';
        S_clu.clusterNotes{iClu2} = '';
    catch
    end
end % function
