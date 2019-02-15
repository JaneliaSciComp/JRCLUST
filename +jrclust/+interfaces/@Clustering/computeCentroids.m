function computeCentroids(obj, updateMe)
    %COMPUTEPOSITIONS determine cluster position from spike position
    if nargin < 2 || isempty(obj.clusterCentroids)
        updateMe = [];
    end

    if isempty(updateMe)
        centroids = zeros(obj.nClusters, 2);
        clusters_ = 1:obj.nClusters;
    else % selective update
        centroids = obj.clusterCentroids;
        clusters_ = updateMe(:)';
    end

    featureSites = 1:obj.hCfg.nSitesEvt;

    for iCluster = clusters_
        [clusterSpikes, neighbors] = subsampleCenteredSpikes(obj, iCluster);
        if isempty(clusterSpikes)
            continue;
        end

        neighbors = neighbors(1:end-obj.hCfg.nSitesExcl);
        featureWeights = squeeze(obj.spikeFeatures(featureSites, 1, clusterSpikes));
        neighborLoc = single(obj.hCfg.siteLoc(neighbors, :)); % position on probe

        centroids(iCluster, 1) = median(getWeightedLoc(featureWeights, neighborLoc(:, 1)));
        centroids(iCluster, 2) = median(getWeightedLoc(featureWeights, neighborLoc(:, 2)));
    end

    obj.clusterCentroids = centroids;
end

%% LOCAL FUNCTIONS
function [clusterSpikes, neighbors] = subsampleCenteredSpikes(hClust, iCluster)
    %SUBSAMPLECENTEREDSPIKES Subsample spikes from the requested cluster centered at the center site and mid-time range (drift)
    nSamplesMax = 1000;

    % subselect based on the center site
    clusterSpikes = hClust.spikesByCluster{iCluster};
    if isempty(clusterSpikes)
        neighbors = [];
        return;
    end

    iClusterSite = hClust.clusterSites(iCluster);
    centeredSpikes = (hClust.spikeSites(clusterSpikes) == iClusterSite);
    neighbors = hClust.hCfg.siteNeighbors(:, iClusterSite);
    clusterSpikes = clusterSpikes(centeredSpikes);
    if isempty(clusterSpikes)
        return;
    end

    clusterSpikes = jrclust.utils.subsample(clusterSpikes, nSamplesMax);
end

function weightedLoc = getWeightedLoc(featureWeights, positions)
    %GETWEIGHTEDLOC Get feature-weighted location of clusters from positions
    if isrow(positions)
        positions = positions';
    end

    featureWeights = featureWeights.^2;
    sos = sum(featureWeights);

    weightedLoc = sum(bsxfun(@times, featureWeights, positions))./sos;
end