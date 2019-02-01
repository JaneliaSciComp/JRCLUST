function spikeData = extractFeatures(spikeData, hCfg)
    %EXTRACTFEATURES Extract spike waveforms and build a spike table
    samplesRaw = spikeData.samplesRaw;
    samplesFilt = spikeData.samplesFilt;
    siteThresh = spikeData.siteThresh;
    spikeTimes = spikeData.spikeTimes;
    spikeSites = spikeData.spikeSites;
    nPadPre = spikeData.nPadPre;

    if hCfg.verbose
        fprintf('\tExtracting features');
        tf = tic;
    end

    spikeSites_ = jrclust.utils.tryGpuArray(spikeSites);

    % spikesRaw, spikesFilt are nSamples x nSites x nSpikes
    [windowsRaw, windowsFilt, spikeTimes] = samplesToWindows(samplesRaw, samplesFilt, spikeSites_, spikeTimes, hCfg);
    if hCfg.verbose
        fprintf('.');
    end

    % get secondary/tertiary peaks and extract windows from them to compute features
    if hCfg.nPeaksFeatures >= 2
        [spikeSites2, spikeSites3] = findSecondaryPeaks(windowsFilt, spikeSites_, hCfg);
        windowsFilt2 = samplesToWindows2(samplesFilt, spikeSites2, spikeTimes, hCfg);

        if hCfg.nPeaksFeatures == 3
            windowsFilt3 = samplesToWindows2(samplesFilt, spikeSites3, spikeTimes, hCfg);
        end
    else
        [spikeSites2, windowsFilt2] = deal([]);
    end

    % Cancel overlap
    if hCfg.getOr('fCancel_overlap', 0)
        try
            [windowsFilt, windowsFilt2] = cancel_overlap_spk_(windowsFilt, windowsFilt2, spikeTimes, spikeSites, spikeSites2, siteThresh, hCfg);
        catch ME
            warning('Cancel overlap failure: %s', ME.message);
        end
    end

    if hCfg.nPeaksFeatures == 1
        features1 = jrclust.features.computeFeatures(windowsFilt, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        spikeFeatures = permute(features1, [1, 3, 2]); % nSites x nFeatures x nSpikes
        centerSites_ = spikeSites_(:);
    elseif hCfg.nPeaksFeatures == 2
        features1 = jrclust.features.computeFeatures(windowsFilt, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        features2 = jrclust.features.computeFeatures(windowsFilt2, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        spikeFeatures = permute(cat(3, features1, features2), [1, 3, 2]); % nSites x nFeatures x nSpikes
        centerSites_ = [spikeSites_(:), spikeSites2(:)]; % nSpikes x nFeatures
    else % hCfg.nPeaksFeatures == 3
        features1 = jrclust.features.computeFeatures(windowsFilt, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        features2 = jrclust.features.computeFeatures(windowsFilt2, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        features3 = jrclust.features.computeFeatures(windowsFilt3, hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        spikeFeatures = permute(cat(3, features1, features2, features3), [1, 3, 2]); % nSites x nFeatures x nSpikes
        centerSites_ = [spikeSites_(:), spikeSites2(:), spikeSites3(:)]; % nSpikes x nFeatures
    end

    if nPadPre > 0
        spikeTimes = spikeTimes - nPadPre;
    end

    spikeData.spikeTimes = jrclust.utils.tryGather(spikeTimes);
    spikeData.spikesRaw = jrclust.utils.tryGather(windowsRaw);
    spikeData.spikesFilt = jrclust.utils.tryGather(windowsFilt);
    spikeData.spikeFeatures = jrclust.utils.tryGather(spikeFeatures);
    spikeData.centerSites = jrclust.utils.tryGather(centerSites_);

    if hCfg.verbose
        fprintf('done (%0.2f s)\n', toc(tf));
    end
end

%% LOCAL FUNCTIONS
function [spikeSites2, spikeSites3] = findSecondaryPeaks(spikeWindows, spikeSites, hCfg)
    %FINDSECONDARYPEAKS Find secondary peaks
    %   spikeWindows: nSamples x nSites x nSpikes
    evtNeighbors = 2:hCfg.nSitesEvt; % sites within merge radius, excluding ref sites
    siteNeighbors2 = hCfg.siteNeighbors(evtNeighbors, spikeSites);

    spikeWindows2 = spikeWindows(:, evtNeighbors, :);
    spikeVpp = squeeze(max(spikeWindows2) - min(spikeWindows2)); % use Vpp to determine second peak site

    [~, peakiestSites] = sort(spikeVpp, 'descend');
    spikeSites2 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(1, :), []));
    spikeSites3 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(2, :), []));
end

function spikeWindows = samplesToWindows2(samplesIn, spikeSites, spikeTimes, hCfg)
    %SAMPLESTOWINDOWS2 Get spatiotemporal windows around secondary peaking
    %events as 3D arrays
    nSpikes = numel(spikeSites);
    nSitesEvt = 1 + hCfg.nSiteDir*2; % includes ref sites

    % nSamples x nSites x nSpikes
    spikeWindows = zeros(diff(hCfg.evtWindowSamp) + 1, nSitesEvt, nSpikes, 'like', samplesIn);

    for iSite = 1:hCfg.nSites
        siteSpikes = find(spikeSites == iSite);
        if isempty(siteSpikes)
            continue;
        end

        siteTimes = spikeTimes(siteSpikes); % already sorted by time
        iNeighbors = hCfg.siteNeighbors(:, iSite);
        spikeWindows(:, :, siteSpikes) = permute(jrclust.utils.tryGather(jrclust.utils.extractWindows(samplesIn, hCfg.evtWindowSamp, siteTimes, iNeighbors)), [1,3,2]);
    end
end
