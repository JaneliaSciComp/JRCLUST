%--------------------------------------------------------------------------
function S_clu = updateSimScore(S_clu)
    % update the cluster similarity score if a kilosort run

    S0 = get(0, 'UserData');
    P = S0.P;

    if getOr(P, 'fImportKilosort', 0)
        rez = S0.rez;

        nClusters = S_clu.nClusters;

        spikeTemplates = S_clu.spikeTemplates;
        spikeClusters = S_clu.spikeClusters;

        simScore = zeros(nClusters);
        for iCluster = 1:nClusters
            spikeClusterIndices = (spikeClusters == iCluster); % spike indices for this cluster
            iClusterTemplates = unique(spikeTemplates(spikeClusterIndices)); % unique template indices for spikes in this cluster

            % compute cluster sim score, Phy style
            sims = max(rez.simScore(iClusterTemplates, :), [], 1);

            for jCluster=iCluster:nClusters
                jClusterTemplates = unique(spikeTemplates(spikeClusters == jCluster));
                simScore(iCluster, jCluster) = max(sims(jClusterTemplates));
                simScore(jCluster, iCluster) = simScore(iCluster, jCluster);
            end
        end

        S_clu.simScore = simScore;
    end
end
