function success = assignClusterCenters(obj)
%ASSIGNCLUSTERCENTERS Summary of this function goes here
%   Detailed explanation goes here
obj.hCfg.updateLog('assignCC', 'Recomputing cluster centers', 1, 0);

success = 0;
clusterCenters = zeros(obj.nClusters, 1);

for iCluster = 1:obj.nClusters
    iSpikes = obj.spikesByCluster{iCluster};

    iRho = obj.spikeRho(iSpikes);
    iDelta = obj.spikeDelta(iSpikes);

    idx = bestClusterCenter(iRho, iDelta);
    clusterCenters(iCluster) = iSpikes(idx);
end

try
    obj.clusterCenters = clusterCenters;
catch ME
    warning(ME.message);
    success = 0;
end

obj.hCfg.updateLog('assignCC', 'Finished recomputing cluster centers', 0, 1);
end

% LOCAL FUNCTIONS
function idx = bestClusterCenter(rho, delta)
%BESTCLUSTERCENTER Find the point with the largest rho value.
%   Ties are broken by the largest delta value.
[maxRho, idx] = max(rho);
maxMask = (rho == maxRho);

if sum(maxMask) > 1
    maxMask = find(maxMask);
    delta = delta(maxMask);
    [~, maxDelta] = max(delta);

    idx = maxMask(maxDelta);
end
end %fun
