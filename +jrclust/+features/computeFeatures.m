function [features1, features2, features3, spikeWindows] = computeFeatures(spikeWindows, hCfg, nSites_spk, spikeSites2)
    %COMPUTEFEATURES Compute features for spikes
    if nargin < 3
        nSites_spk = [];
    end
    if nargin < 4
        spikeSites2 = [];
    end

    features2 = single([]);
    features3 = single([]);
    spikeWindows = single(permute(spikeWindows, [1, 3, 2])); % nSamples x nSpikes x nSites

    if get_set_(hCfg, 'fRealign_spk', 0) ~= 1
        spikeWindows = spkwav_car_(spikeWindows, hCfg, nSites_spk, spikeSites2);
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
    else % pca
        [features1, features2, features3] = jrclust.features.spikePCA(spikeWindows, hCfg);
    end

    if nargout == 1
        if hCfg.nPcPerChan == 2
            features1 = cat(1, features1, features2);
        elseif hCfg.nPcPerChan == 3
            features1 = cat(1, features1, features2, features3);
        end
    end
end