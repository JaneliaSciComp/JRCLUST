function nRemoved = rmRefracSpikes(obj, iCluster)
    %RMREFRACSPIKES remove refractory spikes
    if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

    if nargin == 1 % recurse for each cluster
        nRemoved = 0;
        for iCluster_ = 1:obj.nClusters
            nRemoved_ = obj.rmRefracSpikes(iCluster_);
            nRemoved = nRemoved + nRemoved_;
        end

        return;
    else
        nRemoved = 0;
        nSkipRefrac = 4;

        try
            clusterSpikes_ = obj.spikesByCluster{iCluster};
        catch
            clusterSpikes_ = find(obj.spikeClusters == iCluster);
        end

        if isempty(clusterSpikes_)
            return;
        end

        clusterTimes_ = obj.spikeTimes(clusterSpikes_);

        % removal loop
        keepMe = true(size(clusterTimes_));
        while (1)
            iKeepMe = find(keepMe);

            inRefrac = find(diff(clusterTimes_(keepMe)) < obj.hCfg.refracIntSamp) + 1;
            if isempty(inRefrac)
                break;
            end

            keepMe(iKeepMe(inRefrac(1:nSkipRefrac:end))) = 0;
        end

        nRemoved = sum(~keepMe);
        obj.spikeClusters(clusterSpikes_(~keepMe)) = 0; % assign to noise cluster

        obj.spikesByCluster{iCluster} = clusterSpikes_(keepMe);
        obj.unitCount(iCluster) = sum(keepMe);
    end
end