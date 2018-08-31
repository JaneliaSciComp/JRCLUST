%--------------------------------------------------------------------------
% 171201 JJJ: Unique sites handling for diagonal plotting
function spikeWaveforms = getSpikeWaveformsSites(spikes, sitesOfInterest, S0, fWav_raw_show)
    % get (raw or filtered) spike waveforms on sites of interest
    % return value: nSamples x nSitesOfInterest x nSpikes tensor

    if nargin < 3
        S0 = [];
    end
    if isempty(S0)
        S0 = get(0, 'UserData');
    end
    if nargin < 4
        fWav_raw_show = getOr(S0.P, 'fWav_raw_show', 0);
    end

    % if duplicate sites, operate just on the unique sites
    [uniqueSites, ~, iUniqueSites] = unique(sitesOfInterest);
    if numel(uniqueSites) ~= numel(sitesOfInterest)
        spikeWaveformsUnique = getSpikeWaveformsSites(spikes, uniqueSites, S0, fWav_raw_show);
        spikeWaveforms = spikeWaveformsUnique(:, iUniqueSites, :);
        return;
    end

    P = S0.P;

    % load raw or filtered traces from disk
    traces = getSpikeWaveforms(P, fWav_raw_show);

    nSamples = size(traces, 1);
    nSites = numel(sitesOfInterest);
    nSpikes = numel(spikes);

    spikeSites = S0.spikeSites(spikes);
    spikeSitesUnique = unique(spikeSites);
    spikeWaveforms = zeros([nSamples, nSites, nSpikes], 'like', traces);

    for iSite = 1:numel(spikeSitesUnique)
        thisSite = spikeSitesUnique(iSite);
        iThisSite = find(spikeSites == thisSite);
        neighborsThisSite = P.miSites(:, thisSite);
        [Lia, Locb] = ismember(neighborsThisSite, sitesOfInterest);
        spikeWaveforms(:, Locb(Lia), iThisSite) = traces(:, Lia, spikes(iThisSite));
    end
end % function
