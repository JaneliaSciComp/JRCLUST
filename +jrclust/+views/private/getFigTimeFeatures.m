function [dispFeatures, spikeTimesSecs, YLabel, dispSpikes] = getFigTimeFeatures(hClust, iSite, iCluster)
    %GETFIGTIMEFEATURES Compute features to display in FigTime on iSite
    if nargin < 3 % get features off of background spikes
        iCluster = [];
    end

    hCfg = hClust.hCfg;
    [dispFeatures, dispSpikes] = getClusterFeaturesSite(hClust, iSite, iCluster);
    spikeTimesSecs = double(hClust.spikeTimes(dispSpikes))/hCfg.sampleRate;

    if strcmp(hCfg.dispFeature, 'vpp')
        YLabel = sprintf('Site %d (\\mu Vpp)', iSite);
    else
        YLabel = sprintf('Site %d (%s)', iSite, hCfg.dispFeature);
    end
end

%% LOCAL FUNCTIONS
function [sampledFeatures, sampledSpikes] = getClusterFeaturesSite(hClust, iSite, iCluster)
    %GETCLUSTERFEATURESSITE Get display features for a cluster or
    %background spikes on iSite
    MAX_SAMPLE = 10000; % max points to display

    hCfg = hClust.hCfg;
    spikeSites = hClust.spikeSites;

    if isempty(iCluster) % select spikes based on sites
        nSites = 1 + round(hCfg.nSiteDir);
        neighbors = hCfg.siteNeighbors(1:nSites, iSite);

        if isempty(hClust.spikesBySite)
            sampledSpikes = find(ismember(spikeSites, neighbors));
        else
            sampledSpikes = jrclust.utils.neCell2mat(hClust.spikesBySite(neighbors)');            
        end

        sampledSpikes = jrclust.utils.subsample(sampledSpikes, MAX_SAMPLE);
    else % get all sites from the cluster
        sampledSpikes = hClust.spikesByCluster{iCluster};
    end

    if strcmp(hCfg.dispFeature, 'vpp')
        sampledWaveforms = squeeze(hClust.getSpikeWindows(sampledSpikes, iSite, 0, 1)); % use voltages
        sampledFeatures = max(sampledWaveforms) - min(sampledWaveforms);
    elseif strcmp(hCfg.dispFeature, 'cov')
        sampledFeatures = getSpikeCov(hClust, sampledSpikes, iSite);
    elseif strcmp(hCfg.dispFeature, 'pca') || (strcmp(hCfg.dispFeature, 'ppca') && isempty(iCluster)) % TODO: need a better mech for bg spikes
        sampledWindows = permute(hClust.getSpikeWindows(sampledSpikes, iSite, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
        prVecs1 = jrclust.features.getPVSpikes(sampledWindows);
        sampledFeatures = jrclust.features.pcProjectSpikes(sampledWindows, prVecs1);
    elseif strcmp(hCfg.dispFeature, 'ppca')
        sampledWindows = permute(hClust.getSpikeWindows(sampledSpikes, iSite, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
        prVecs1 = jrclust.features.getPVClusters(hClust, iSite, iCluster);
        sampledFeatures = jrclust.features.pcProjectSpikes(sampledWindows, prVecs1);
    elseif strcmp(hCfg.dispFeature, 'template')
        sampledFeatures = hClust.templateFeaturesBySpike(sampledSpikes, iSite);
    else
        error('not implemented yet');
    end

    sampledFeatures = squeeze(abs(sampledFeatures));
end

