function success = commit(obj, spikeClusters, metadata, msg)
    %COMMIT Commit a modification of clustering to history log
    success = commit@jrclust.interfaces.Clustering(obj, spikeClusters, metadata, msg);

    if success && isfield(metadata, 'unitCount') && any(isnan(metadata.unitCount))
        flagged = find(isnan(metadata.unitCount));

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

