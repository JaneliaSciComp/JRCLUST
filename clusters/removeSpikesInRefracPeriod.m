%--------------------------------------------------------------------------
function [S_clu, nRemoved] = removeSpikesInRefracPeriod(S_clu, P, cluster)
    % remove spikes in refractory period

    spikeTimes = get0_('spikeTimes');

    if nargin == 2
        nClusters = max(S_clu.spikeClusters);
        nRemoved = 0;

        for iCluster = 1:nClusters
            [S_clu, iRemoved] = removeSpikesInRefracPeriod(S_clu, P, iCluster);
            nRemoved = nRemoved + iRemoved;
        end

        return;
    else
        nRemoved = 0;
        nSkip_refrac = getOr(P, 'nSkip_refrac', 4);

        try
            spikes = S_clu.spikesByCluster{cluster};
        catch
            spikes = find(S_clu.spikeClusters == cluster);
        end

        if isempty(spikes)
            return;
        end

        clusterTimes = spikeTimes(spikes);
        nRefrac = ceil(P.spkRefrac_ms * P.sampleRateHz / 1000);

        % removal loop
        spikesToKeep = true(size(spikes));
        while 1
            iSpikesToKeep = find(spikesToKeep);

            % find spikes occurring too close together in time
            violations = find(diff(clusterTimes(iSpikesToKeep)) < nRefrac) + 1;
            if isempty(violations)
                break;
            end

            spikesToKeep(iSpikesToKeep(violations(1:nSkip_refrac:end))) = 0;
        end

        nRemoved = sum(~spikesToKeep);
        nSpikes = numel(spikesToKeep);

        S_clu.spikeClusters(spikes(~spikesToKeep)) = 0;

        S_clu.spikesByCluster{cluster} = spikes(spikesToKeep);
        S_clu.nSpikesPerCluster(cluster) = sum(spikesToKeep);
    end

    if get_(P, 'fVerbose')
        fprintf('Removed %d/%d (%0.1f%%) duplicate spikes from cluster %d\n', nRemoved, nSpikes, nRemoved/nSpikes*100, cluster);
    end
end % function
