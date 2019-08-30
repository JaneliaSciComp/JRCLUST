function flag = recover(obj, ask)
    %RECOVER Inconsistent data recovery, with a hammer
    if nargin < 2
        ask = 0;
    end

    flag = recover@jrclust.interfaces.Clustering(obj, ask);
    if flag < 1
        return;
    end

    % update cluster centers
    if flag == 2 || ismember('clusterCenters', obj.inconsistentFields())
        for iCluster = 1:numel(obj.nClusters)
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

        % truncate "extra" centers, if any
        obj.clusterCenters = obj.clusterCenters(1:obj.nClusters);
    end
end

