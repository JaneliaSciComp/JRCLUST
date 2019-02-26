function spikeData = extractFeatures(obj, spikeData)
    %EXTRACTFEATURES Extract spike waveforms and build a spike table
    spikesFilt = jrclust.utils.tryGpuArray(spikeData.spikesFilt, obj.hCfg.useGPU);
    spikesFilt2 = jrclust.utils.tryGpuArray(spikeData.spikesFilt2, obj.hCfg.useGPU);
    spikesFilt3 = jrclust.utils.tryGpuArray(spikeData.spikesFilt3, obj.hCfg.useGPU);

    obj.hCfg.updateLog('extractFeatures', 'Extracting features', 1, 0);

    % spikesRaw, spikesFilt are nSamples x nSites x nSpikes
    %[windowsRaw, windowsFilt, spikeTimes] = obj.samplesToWindows(samplesRaw, samplesFilt, spikeTimes, spikeSites_);

    % get secondary/tertiary peaks and extract windows from them to compute features
%     if obj.hCfg.nPeaksFeatures >= 2
%         [spikeSites2, spikeSites3] = obj.findSecondaryPeaks(windowsFilt, spikeSites_);
%         windowsFilt2 = obj.samplesToWindows2(samplesFilt, spikeSites2, spikeTimes);
% 
%         if obj.hCfg.nPeaksFeatures == 3
%             windowsFilt3 = obj.samplesToWindows2(samplesFilt, spikeSites3, spikeTimes);
%         end
%     else
%         [spikeSites2, windowsFilt2] = deal([]);
%     end

    if obj.hCfg.nPeaksFeatures == 1
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);

        spikeFeatures = permute(features1, [1, 3, 2]); % nSites x nFeatures x nSpikes
    elseif obj.hCfg.nPeaksFeatures == 2
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        features2 = jrclust.features.computeFeatures(spikesFilt2, obj.hCfg);

        spikeFeatures = permute(cat(3, features1, features2), [1, 3, 2]); % nSites x nFeatures x nSpikes
    else % obj.hCfg.nPeaksFeatures == 3
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        features2 = jrclust.features.computeFeatures(spikesFilt2, obj.hCfg);
        features3 = jrclust.features.computeFeatures(spikesFilt3, obj.hCfg);

        spikeFeatures = permute(cat(3, features1, features2, features3), [1, 3, 2]); % nSites x nFeatures x nSpikes
    end

    [spikesFilt, spikesFilt2, spikesFilt3, spikeFeatures] = jrclust.utils.tryGather(spikesFilt, spikesFilt2, spikesFilt3, spikeFeatures); %#ok<ASGLU>
    spikeData.spikeFeatures = spikeFeatures;
    obj.hCfg.updateLog('extractFeatures', 'Finished extracting features', 0, 1);
end
