function S_clu_position_(hClust, updateMe)
    % determine cluster position from spike position
    % 6/27/17 JJJ: multiple features supported (single dimension such as energy and Vpp)
    if nargin < 2 || isempty(hClust.centroids)
        updateMe = [];
    end

    if isempty(updateMe)
        centroids_ = zeros(hClust.nClusters, 2);
        clusters_ = 1:hClust.nClusters;
    else % selective update
        centroids_ = hClust.clusterCentroids;
        clusters_ = updateMe(:)';
    end

    featureSites = 1:hClust.hCfg.nSitesEvt);

    for iCluster = clusters_
        [clusterSpikes, neighbors] = subsampleCenteredSpikes(hClust, iCluster);
        if isempty(clusterSpikes)
            continue;
        end

        neighbors = neighbors(1:end-hClust.hCfg.nSitesExcl);
        featureWeights = squeeze(hClust.spikeFeatures(featureSites, 1, clusterSpikes));
        neighborLoc = single(hClust.hCfg.siteLoc(neighbors, :)); % position on probe

        centroids_(1, iCluster) = median(getWeightedLoc(featureWeights, neighborLoc(:, 1)));
        centroids_(2, iCluster) = median(getWeightedLoc(featureWeights, neighborLoc(:, 2)));
    end

    hClust.clusterCentroids = centroids_;
end