%function tnWav_spk1 = jrclust.utils.getSampledWindows(sampledSpikes, iSite, S0, fWav_raw_show)
function sampledWindows = getSampledWindows(hClust, sampledSpikes, iSite, showRaw)
    hCfg = hClust.hCfg;
    spikeSites = hClust.spikeSites;

    if nargin < 4
        showRaw = hCfg.getOr('showRaw', false);
    end

    % don't repeat ourselves if we have duplicate sites
    [uniqueSites, ~, indMap] = unique(iSite);
    if numel(uniqueSites) ~= numel(iSite)
        swUnique = getSampledWindows(sampledSpikes, uniqueSites, S0, showRaw);
        sampledWindows = swUnique(:, indMap, :);
        return;
    end

    if showRaw
        spikeWindows = hClust.spikesRaw;
    else
        spikeWindows = hClust.spikesFilt;
    end
    nSamples = size(spikeWindows, 1);
    nSpikes = numel(sampledSpikes);

    sampledSites = spikeSites(sampledSpikes);
    uniqueSites = unique(sampledSites);

    sampledWindows = zeros([nSamples, numel(iSite), nSpikes], 'like', spikeWindows);
    for jSite = 1:numel(uniqueSites)
        site = uniqueSites(jSite); % center sites group
        neighbors = hCfg.siteNeighbors(:, site);

        onSite = find(sampledSites == site);
        [isNeigh, neighLocs] = ismember(neighbors, iSite);
        sampledWindows(:, neighLocs(isNeigh), onSite) = spikeWindows(:, isNeigh, sampledSpikes(onSite));
    end
end
