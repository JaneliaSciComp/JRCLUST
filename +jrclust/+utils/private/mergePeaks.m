function [spikeTimes, spikeAmps, spikeSites] = mergePeaks(spikesBySite, ampsBySite, hCfg)
    %MERGEPEAKS Merge duplicate peak events
    nSites = numel(spikesBySite);
    spikeTimes = jrclust.utils.neCell2mat(spikesBySite);
    spikeAmps = jrclust.utils.neCell2mat(ampsBySite);
    spikeSites = jrclust.utils.neCell2mat(cellfun(@(vi, i) repmat(i, size(vi)), spikesBySite, num2cell((1:nSites)'), 'UniformOutput', 0));

    [spikeTimes, argsort] = sort(spikeTimes);
    spikeAmps = spikeAmps(argsort);
    spikeSites = int32(spikeSites(argsort));
    spikeTimes = int32(spikeTimes);

    [mergedTimes, mergedAmps, mergedSites] = deal(cell(nSites,1));

    try
        % avoid sending the entire hCfg object out to workers
        cfgSub = struct('refracIntSamp', hCfg.refracIntSamp, ...
                        'siteLoc', obj.hCfg.siteLoc, ...
                        'evtDetectRad', obj.hCfg.evtDetectRad);
        parfor iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, cfgSub);
            catch % don't try to display an error here
            end
        end
    catch % parfor failure
        for iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch ME
                warning('failed to merge spikes on site %d: %s', iSite, ME.message);
            end
        end
    end

    % merge parfor output and sort
    spikeTimes = jrclust.utils.neCell2mat(mergedTimes);
    spikeAmps = jrclust.utils.neCell2mat(mergedAmps);
    spikeSites = jrclust.utils.neCell2mat(mergedSites);

    [spikeTimes, argsort] = sort(spikeTimes); % sort by time
    spikeAmps = jrclust.utils.tryGather(spikeAmps(argsort));
    spikeSites = spikeSites(argsort);
end

%% LOCAL FUNCTIONS
function [timesOut, ampsOut, sitesOut] = mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg)
    %MERGESPIKESSITE Merge spikes in the refractory period
    nLims = int32(abs(hCfg.refracIntSamp));

    % find neighboring spikes
    nearbySites = jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtDetectRad); % includes iSite
    spikesBySite = arrayfun(@(jSite) find(spikeSites == jSite), nearbySites, 'UniformOutput', 0);
    timesBySite = arrayfun(@(jSite) spikeTimes(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);
    ampsBySite = arrayfun(@(jSite) spikeAmps(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);

    iiSite = (nearbySites == iSite);
    iSpikes = spikesBySite{iiSite};
    iTimes = timesBySite{iiSite};
    iAmps = ampsBySite{iiSite};

    % search over peaks on neighboring sites and in refractory period to
    % see which peaks on this site to keep
    keepMe = true(size(iSpikes));
    for jjSite = 1:numel(spikesBySite)
        jSite = nearbySites(jjSite);

        jSpikes = spikesBySite{jjSite};
        jTimes = timesBySite{jjSite};
        jAmps = ampsBySite{jjSite};

        if iSite == jSite
            delays = [-nLims:-1, 1:nLims]; % skip 0 delay
        else
            delays = -nLims:nLims;
        end

        for iiDelay = 1:numel(delays)
            iDelay = delays(iiDelay);

            [jLocs, iLocs] = ismember(jTimes, iTimes + iDelay);
            if ~any(jLocs)
                continue;
            end

            % jLocs: index into j{Spikes,Times,Amps}
            % iLocs: (1st) corresponding index into i{Spikes,Times,Amps}/keepMe
            jLocs = find(jLocs);
            iLocs = iLocs(jLocs);

            % drop spikes on iSite where spikes on jSite have larger
            % magnitudes
            nearbyLarger = abs(jAmps(jLocs)) > abs(iAmps(iLocs));
            keepMe(iLocs(nearbyLarger)) = 0;

            % flag equal-amplitude nearby spikes
            ampsEqual = (jAmps(jLocs) == iAmps(iLocs));
            if any(ampsEqual)
                if iDelay < 0 % drop spikes on iSite occurring later than those on jSite
                    keepMe(iLocs(ampsEqual)) = 0;
                elseif iDelay == 0 && jSite < iSite % drop spikes on iSite if jSite is of lower index
                    keepMe(iLocs(ampsEqual)) = 0;
                end
            end
        end
    end

    % keep the peak spikes only
    timesOut = iTimes(keepMe);
    ampsOut = iAmps(keepMe);
    sitesOut = repmat(int32(iSite), size(timesOut));
end
