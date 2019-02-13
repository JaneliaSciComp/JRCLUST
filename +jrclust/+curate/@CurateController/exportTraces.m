function exportTraces(obj)
    %EXPORTTRACES Export traces from selected cluster to workspace
    spikesFilt = obj.hClust.spikesFilt;
    spikesRaw = obj.hClust.spikesRaw;
    iCluster = obj.selected(1);
    clusterSpikes = obj.hClust.spikesByCluster{iCluster};
    nSpikes = numel(clusterSpikes);

    shapeFilt = [size(spikesFilt, 1), obj.hCfg.nSites, nSpikes];
    shapeRaw = [size(spikesRaw, 1), obj.hCfg.nSites, nSpikes];

    iSpikesFilt = zeros(shapeFilt, 'like', spikesFilt);
    iSpikesRaw = zeros(shapeRaw, 'like', spikesRaw);
    spikeNeighbors = obj.hCfg.siteNeighbors(:, obj.hClust.spikeSites);

    for jSpike = 1:nSpikes
        iSpikesFilt(:, spikeNeighbors(:, jSpike), jSpike) = spikesFilt(:, :, clusterSpikes(jSpike));
        iSpikesRaw(:, spikeNeighbors(:, jSpike), jSpike) = spikesRaw(:, :, clusterSpikes(jSpike));
    end

    jrclust.utils.exportToWorkspace(struct(sprintf('spikesFilt%d', iCluster), iSpikesFilt, ...
                                           sprintf('spikesRaw%d', iCluster), iSpikesRaw), obj.hCfg.verbose);
end