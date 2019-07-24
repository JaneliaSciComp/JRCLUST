function success = revert(obj, revertTo)
    %REVERT Summary of this function goes here
    success = revert@jrclust.interfaces.Clustering(obj, revertTo);

    if success
        good = find(obj.clusterCenters > 0);
        % find both empty cluster centers and centers that don't belong to
        % this cluster
        flagged = [find(obj.clusterCenters == 0); good(obj.spikeClusters(obj.clusterCenters(good)) ~= good)];

        % update cluster centers
        for i = 1:numel(flagged)
            iCluster = flagged(i);
            iSpikes = obj.spikesByCluster{iCluster};
            iRho = obj.spikeRho(iSpikes);
            [maxRho, newCenter] = max(iRho);
            maxMask = (iRho == maxRho);

            % tie in spike rho, use delta as a tie breaker
            if sum(maxMask) > 1
                iDelta = obj.spikeDelta(iSpikes);
                maxMask = find(maxMask);

                [~, maxDelta] = max(iDelta(maxMask));
                newCenter = maxMask(maxDelta);
            end

            obj.clusterCenters(iCluster) = iSpikes(newCenter);
        end
    end
end

