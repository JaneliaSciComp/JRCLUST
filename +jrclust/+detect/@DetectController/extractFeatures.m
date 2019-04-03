function spikeData = extractFeatures(obj, spikeData)
    %EXTRACTFEATURES Extract features from filtered spike traces
    if obj.hCfg.useGPU
        S = gpuDevice();
    end

    obj.hCfg.updateLog('extractFeatures', 'Extracting features', 1, 0);

    if obj.hCfg.nPeaksFeatures >= 1
        % try to compute in GPU iff we have the memory
        if obj.hCfg.useGPU && floor(log10(numel(spikeData.spikesFilt)*jrclust.utils.typeBytes('single'))) < floor(log10(S(1).AvailableMemory))
            spikesFilt = jrclust.utils.tryGpuArray(spikeData.spikesFilt, obj.hCfg.useGPU);
        else
            spikesFilt = spikeData.spikesFilt;
        end

        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        [spikesFilt, features1] = jrclust.utils.tryGather(spikesFilt, features1); %#ok<*ASGLU>

        %spikeFeatures = permute(features1, [1, 3, 2]); % nSites x nFeatures x nSpikes
        spikeFeatures = features1;
    end

    if obj.hCfg.nPeaksFeatures >= 2
        % try to compute in GPU iff we have the memory
        if obj.hCfg.useGPU && floor(log10(numel(spikeData.spikesFilt2)*jrclust.utils.typeBytes('single'))) < floor(log10(S(1).AvailableMemory))
            spikesFilt2 = jrclust.utils.tryGpuArray(spikeData.spikesFilt2, obj.hCfg.useGPU);
        else
            spikesFilt2 = spikeData.spikesFilt2;
        end

        features2 = jrclust.features.computeFeatures(spikesFilt2, obj.hCfg);
        [spikesFilt2, features2] = jrclust.utils.tryGather(spikesFilt2, features2);

        %spikeFeatures = permute(cat(3, features1, features2), [1, 3, 2]); % nSites x nFeatures x nSpikes
        spikeFeatures = cat(3, spikeFeatures, features2);
    end

    if obj.hCfg.nPeaksFeatures == 3
        % try to compute in GPU iff we have the memory
        if obj.hCfg.useGPU && floor(log10(numel(spikeData.spikesFilt3)*jrclust.utils.typeBytes('single'))) < floor(log10(S(1).AvailableMemory))
            spikesFilt3 = jrclust.utils.tryGpuArray(spikeData.spikesFilt3, obj.hCfg.useGPU);
        else
            spikesFilt3 = spikeData.spikesFilt3;
        end

        features3 = jrclust.features.computeFeatures(spikesFilt3, obj.hCfg);
        [spikesFilt3, features3] = jrclust.utils.tryGather(spikesFilt3, features3);

        %spikeFeatures = permute(cat(3, features1, features2, features3), [1, 3, 2]); % nSites x nFeatures x nSpikes
        spikeFeatures = cat(3, spikeFeatures, features3);
    end

    spikeData.spikeFeatures = permute(spikeFeatures, [1, 3, 2]); % nSites x nFeatures x nSpikes
    obj.hCfg.updateLog('extractFeatures', 'Finished extracting features', 0, 1);
end
