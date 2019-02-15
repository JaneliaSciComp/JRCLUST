function sampledWindows = getSpikeWindows(obj, spikes, sites, useRaw, useVolt)
    %GETSPIKEWINDOWS 
    %   input:  spikes, a vector of spike indices
    %   input:  sites, a vector of site indices
    %   input:  useRaw, bool, a flag use raw waveforms if true (else filtered)
    %   input:  useVolt, bool, a flag to use traces in units of microvolts
    %           if true (else use the digitized values)
    %   output: sampledWindows is nSamples x nSites x nSpikes
    hCfg = obj.hCfg;
    if nargin < 4
        useRaw = hCfg.getOr('showRaw', 0);
    end
    if nargin < 5
        useVolt = (useRaw && ~isempty(obj.spikesRawVolt)) || (~useRaw && ~isempty(obj.spikesFiltVolt));
    end

    % don't repeat ourselves if we have duplicate sites
    [uniqueSites, ~, inds] = unique(sites);
    if numel(uniqueSites) ~= numel(sites)
        swUnique = obj.getSpikeWindows(spikes, uniqueSites, S0, useRaw, useVolt);
        sampledWindows = swUnique(:, inds, :);
        return;
    end

    if useRaw
        if useVolt
            if isempty(obj.spikesRawVolt) % compute on the fly and cache
                obj.spikesRawVolt = jrclust.utils.rawTouV(obj.spikesRaw, obj.hCfg);
            end

            spikeWindows = obj.spikesRawVolt;
        else
            spikeWindows = obj.spikesRaw;
        end
    else
        if useVolt
            if isempty(obj.spikesFiltVolt) % compute on the fly and cache
                obj.spikesFiltVolt = jrclust.utils.filtTouV(obj.spikesFilt, obj.hCfg);
            end

            spikeWindows = obj.spikesFiltVolt;
        else
            spikeWindows = obj.spikesFilt;
        end
    end

    nSamples = size(spikeWindows, 1);
    nSpikes = numel(spikes);

    sampledSites = obj.spikeSites(spikes);
    uniqueSites = unique(sampledSites);

    sampledWindows = zeros([nSamples, numel(sites), nSpikes], 'like', spikeWindows);
    for jSite = 1:numel(uniqueSites)
        site = uniqueSites(jSite); % center sites group
        neighbors = hCfg.siteNeighbors(:, site);

        onSite = find(sampledSites == site);
        [isNeigh, neighLocs] = ismember(neighbors, sites);
        sampledWindows(:, neighLocs(isNeigh), onSite) = spikeWindows(:, isNeigh, spikes(onSite));
    end
end
