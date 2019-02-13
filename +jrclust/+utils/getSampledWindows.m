function sampledWindows = getSampledWindows(hClust, sampledSpikes, sites, showRaw)
    %GETSAMPLEDWINDOWS 
    %   input:  sampledSpikes, a vector of spike indices
    %   input:  sites, a vector of site indices
    %   input:  showRaw, bool, a flag use raw waveforms if true (else filtered)
    %   output: sampledWindows is nSamples x nSites x nSpikes
    hCfg = hClust.hCfg;
    if nargin < 4
        showRaw = hCfg.getOr('showRaw', 0);
    end

    % don't repeat ourselves if we have duplicate sites
    [uniqueSites, ~, inds] = unique(sites);
    if numel(uniqueSites) ~= numel(sites)
        swUnique = getSampledWindows(sampledSpikes, uniqueSites, S0, showRaw);
        sampledWindows = swUnique(:, inds, :);
        return;
    end

    if showRaw
        spikeWindows = hClust.spikesRaw;
    else
        spikeWindows = hClust.spikesFilt;
    end

    nSamples = size(spikeWindows, 1);
    nSpikes = numel(sampledSpikes);

    sampledSites = hClust.spikeSites(sampledSpikes);
    uniqueSites = unique(sampledSites);

    sampledWindows = zeros([nSamples, numel(sites), nSpikes], 'like', spikeWindows);
    for jSite = 1:numel(uniqueSites)
        site = uniqueSites(jSite); % center sites group
        neighbors = hCfg.siteNeighbors(:, site);

        onSite = find(sampledSites == site);
        [isNeigh, neighLocs] = ismember(neighbors, sites);
        sampledWindows(:, neighLocs(isNeigh), onSite) = spikeWindows(:, isNeigh, sampledSpikes(onSite));
    end
end
