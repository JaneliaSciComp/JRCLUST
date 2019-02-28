function [features1, features2, features3, spikeWindows] = computeFeatures(spikeWindows, hCfg, nSitesEvt)
    %COMPUTEFEATURES Compute features for spikes
    if nargin < 3
        nSitesEvt = [];
    end

    [features2, features3] = deal(single([]));

    spikeWindows = single(permute(spikeWindows, [1, 3, 2])); % nSamples x nSpikes x nSites

    if hCfg.getOr('realignTraces', 0) ~= 1
        spikeWindows = jrclust.utils.localCAR(spikeWindows, hCfg, nSitesEvt, []);
    end

    if strcmp(hCfg.clusterFeature, 'cov')
        [features1, features2] = jrclust.features.spikeCov(spikeWindows, hCfg);
    elseif strcmp(hCfg.clusterFeature, 'vpp')
        [~, ~, features1] = jrclust.features.spikeMinMax(spikeWindows);
    elseif strcmp(hCfg.clusterFeature, 'vmin')
        features1 = jrclust.features.spikeMinMax(spikeWindows);
    elseif strcmp(hCfg.clusterFeature, 'vminmax')
        [features1, features2] = jrclust.features.spikeMinMax(spikeWindows);
    elseif strcmp(hCfg.clusterFeature, 'energy')
        features1 = jrclust.features.spikeEnergy(spikeWindows);
    else % pca, gpca
        [features1, features2, features3] = jrclust.features.spikePCA(spikeWindows, hCfg);
    end

    if nargout == 1
        if hCfg.nPCsPerSite == 2
            features1 = cat(1, features1, features2);
        elseif hCfg.nPCsPerSite == 3
            features1 = cat(1, features1, features2, features3);
        end
    end
end