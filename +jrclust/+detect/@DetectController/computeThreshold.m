function siteThresh = computeThreshold(obj, samplesIn)
    %COMPUTETHRESHOLD Compute sitewise threshold for samplesIn
    maxSample = 1e5;

    try
        samplesIn = jrclust.utils.tryGpuArray(samplesIn, obj.hCfg.useGPU);
        siteThresh = jrclust.utils.estimateRMS(samplesIn, maxSample)*obj.hCfg.qqFactor;
        [samplesIn, siteThresh] = jrclust.utils.tryGather(samplesIn, siteThresh); %#ok<ASGLU>
    catch ME
        warning('GPU threshold computation failed: %s (retrying in CPU)', ME.message);
        samplesIn = jrclust.utils.tryGather(samplesIn);
        obj.hCfg.useGPU = 0;
        siteThresh = jrclust.utils.estimateRMS(samplesIn, maxSample)*obj.hCfg.qqFactor;
    end

    siteThresh = int16(siteThresh);
end
