function [spikeTimes, spikeAmps, spikeSites] = mergePeaks(spikesBySite, ampsBySite, hCfg)
    %MERGEPEAKS Merge duplicate peak events
    nSites = numel(spikesBySite);
    spikeTimes = jrclust.utils.neCell2mat(spikesBySite);
    spikeAmps = jrclust.utils.neCell2mat(ampsBySite);
    spikeSites = jrclust.utils.neCell2mat(cellfun(@(vi, i) repmat(i, size(vi)), spikesBySite, num2cell((1:nSites)'), 'UniformOutput', false));

    [spikeTimes, argsort] = sort(spikeTimes);
    spikeAmps = spikeAmps(argsort);
    spikeSites = int32(spikeSites(argsort));
    spikeTimes = int32(spikeTimes);

    [mergedTimes, mergedAmps, mergedSites] = deal(cell(nSites,1));

    try
        parfor iSite = 1:nSites %parfor speedup: 2x %parfor
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch
                disperr_();
            end
        end
    catch % parfor failure
        for iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch
                disperr_();
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

    % find spikes on iSite
    onSite = int32(find(spikeSites == iSite)); % pre-cache
    onSiteTimes = spikeTimes(onSite);
    onSiteAmps = spikeAmps(onSite);

    % find neighboring spikes
    nearbySites = jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtMergeRad); % includes iSite
    onNearSite = int32(find(ismember(spikeSites, nearbySites)));
    neighborTimes = spikeTimes(onNearSite);
    neighborAmps = spikeAmps(onNearSite);
    neighborSites = spikeSites(onNearSite);

    % search over peaks on neighboring sites and in refractory period to
    % see which peaks on this site to keep
    keepMe = true(size(onSite));
    for iDelay = -nLims:nLims
        [isNearby, locNearby] = ismember(neighborTimes, onSiteTimes + iDelay);
        isNearby = find(isNearby);

        if iDelay == 0 % remove self if zero delay
            isNearby(onNearSite(isNearby) == onSite(locNearby(isNearby))) = [];
        end

        % ignore nearby spikes which have smaller ("less negative") magnitudes
        isNearby(neighborAmps(isNearby) > onSiteAmps(locNearby(isNearby))) = [];
        
        % flag equal-amplitude nearby spikes
        ampsEqual = (neighborAmps(isNearby) == onSiteAmps(locNearby(isNearby)));
        if any(ampsEqual)
            if iDelay > 0 % spk1 occurs before spk12, thus keep
                isNearby(ampsEqual) = [];
            elseif iDelay == 0 % keep only if site is lower
                ampsEqual(iSite > neighborSites(isNearby(ampsEqual))) = 0;
                isNearby(ampsEqual) = []; % same site/time/amp
            end
        end

        % discard spikes on this site which have not passed our tests
        keepMe(locNearby(isNearby)) = 0;
    end %for

    % keep the peak spikes only
    timesOut = onSiteTimes(keepMe);
    ampsOut = onSiteAmps(keepMe);
    sitesOut = repmat(int32(iSite), size(timesOut));
end