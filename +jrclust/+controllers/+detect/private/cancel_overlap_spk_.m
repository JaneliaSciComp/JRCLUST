% 12/16/17 JJJ: Find overlapping spikes and set superthreshold sample points to zero in the overlapping region
function [tnWav_spk_out, tnWav_spk2_out] = cancel_overlap_spk_(spikeWindows, spikeWindows2, spikeTimes, spikeSites, spikeSites2, siteThresh, hCfg)
    % Overlap detection. only return one stronger than other
    [spikeTimes, spikeWindows, spikeWindows2] = jrclust.utils.tryGather(spikeTimes, spikeWindows, spikeWindows2);
    [viSpk_ol_spk, vnDelay_ol_spk] = findPotentialOverlaps(spikeTimes, spikeSites, hCfg);
    [tnWav_spk_out, tnWav_spk2_out] = deal(spikeWindows, spikeWindows2);

    % find spike index that are larger and fit and deploy
    viSpk_ol_a = find(viSpk_ol_spk>0); % later occuring
    [viSpk_ol_b, vnDelay_ol_b] = deal(viSpk_ol_spk(viSpk_ol_a), vnDelay_ol_spk(viSpk_ol_a)); % first occuring
    viTime_spk0 = int32(hCfg.evtWindowSamp(1):hCfg.evtWindowSamp(2));
    siteThresh = jrclust.utils.tryGather(-abs(siteThresh(:))');

    % for each pair identify time range where threshold crossing occurs and set to zero
    % correct only first occuring (b)
    siteNeighbors = hCfg.siteNeighbors;
    nSpk_ol = numel(viSpk_ol_a);
    nSpk = size(spikeWindows,2);

    for iSpk_ol = 1:nSpk_ol
        [iSpk_b, nDelay_b] = deal(viSpk_ol_b(iSpk_ol), vnDelay_ol_b(iSpk_ol));
        viSite_b = siteNeighbors(:,spikeSites(iSpk_b));
        mnWav_b = tnWav_spk_out(nDelay_b+1:end,:,iSpk_b);
        mlWav_b = bsxfun(@le, mnWav_b, siteThresh(viSite_b));
        mnWav_b(mlWav_b) = 0;
        tnWav_spk_out(nDelay_b+1:end,:,iSpk_b) = mnWav_b;

        if ~isempty(spikeWindows2)
            viSite_b = siteNeighbors(:,spikeSites2(iSpk_b));
            mnWav_b = tnWav_spk2_out(nDelay_b+1:end,:,iSpk_b);
            mlWav_b = bsxfun(@le, mnWav_b, siteThresh(viSite_b));
            mnWav_b(mlWav_b) = 0;
            tnWav_spk2_out(nDelay_b+1:end,:,iSpk_b) = mnWav_b;
        end
    end

    % set no overthreshold zone based on the delay, set it to half. only set superthreshold spikes to zero
end

%% LOCAL FUNCTIONS
function [viSpk_ol_spk, vnDelay_ol_spk] = findPotentialOverlaps(spikeTimes, spikeSites, hCfg)
    %FINDPOTENTIALOVERLAPS
    nSites = max(spikeSites);
    spikesBySite = arrayfun(@(iSite) int32(find(spikeSites == iSite)), 1:nSites, 'UniformOutput', 0);
    spikeTimes = jrclust.utils.tryGather(spikeTimes);
    [viSpk_ol_spk, vnDelay_ol_spk] = deal(zeros(size(spikeSites), 'int32'));

    for iSite = 1:nSites
        spikesOnSite = spikesBySite{iSite}; % spikes occurring on this site
        if isempty(spikesOnSite)
            continue;
        end

        % get sites which are *near* iSite but not iSite
        nearbySites = setdiff(jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtDetectRad), iSite);
        nearbySpikes = jrclust.utils.neCell2mat(spikesBySite(nearbySites));

        nSpikesOnSite = numel(spikesOnSite);
        neighborSpikes = [spikesOnSite(:); nearbySpikes(:)]; % spikes in a neighborhood of iSite (includes iSite)

        onSiteTimes = spikeTimes(spikesOnSite);
        neighborTimes = spikeTimes(neighborSpikes);

        % search over the event window for spikes occurring too close in
        % space and time
        for iDelay = 0:diff(hCfg.evtWindowSamp)
            [isNearby, locNearby] = ismember(neighborTimes, onSiteTimes + iDelay);
            if iDelay == 0
                isNearby(1:nSpikesOnSite) = 0; % exclude same site comparison
            end

            isNearby = find(isNearby);
            if isempty(isNearby)
                continue;
            end

            viSpk1_ = spikesOnSite(locNearby(isNearby));
            viSpk12_ = neighborSpikes(isNearby);

            viSpk_ol_spk(viSpk12_) = viSpk1_;
            vnDelay_ol_spk(viSpk12_) = iDelay;
        end
    end
end
