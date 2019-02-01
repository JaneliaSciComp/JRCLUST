%function [mrVpp1, mrVpp2] = getSpikeCov(sampledSpikes, iSite)
function [mrVpp1, mrVpp2] = getSpikeCov(hClust, sampledSpikes, iSite)
    %GETSPIKECOV Compute covariance feature for a subset of spikes
    hCfg = hClust.hCfg;
    spikeSites = hClust.spikeSites;
    spikesFilt = hClust.spikesFilt;

    nSpikes = numel(sampledSpikes);
    sampledSites = spikeSites(sampledSpikes);

    % nSamples x nSpikes x nSites
    sampledWindows = single(permute(spikesFilt(:, :, sampledSpikes), [1, 3, 2]));
    sampledWindows = jrclust.utils.tryGpuArray(sampledWindows, hCfg.useGPU);

    [mrVpp1_, mrVpp2_] = jrclust.features.spikeCov(sampledWindows, hCfg);
    mrVpp1_ = jrclust.utils.tryGather(abs(mrVpp1_));
    mrVpp2_ = jrclust.utils.tryGather(abs(mrVpp2_));

    % re-project to common basis
    uniqueSites = unique(sampledSites);
    [mrVpp1, mrVpp2] = deal(zeros([numel(iSite), nSpikes], 'like', mrVpp1_));

    for jSite = 1:numel(uniqueSites)
        site = uniqueSites(jSite);
        neighbors = hCfg.siteNeighbors(:, site);
        onSite = find(sampledSites == site);

        [isNeigh, neighLocs] = ismember(neighbors, iSite);
        mrVpp1(neighLocs(isNeigh), onSite) = mrVpp1_(isNeigh, onSite);
        mrVpp2(neighLocs(isNeigh), onSite) = mrVpp2_(isNeigh, onSite);
    end
end
