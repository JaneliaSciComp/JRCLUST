function siteThresh = computeThreshold(obj, samplesIn)
    %COMPUTETHRESHOLD Compute sitewise threshold for samplesIn
    try
        siteThresh = jrclust.utils.estimateRMS(samplesIn, 1e5)*obj.hCfg.qqFactor;
        siteThresh = int16(jrclust.utils.tryGather(siteThresh));
    catch ME
        warning('GPU threshold computation failed: %s (retrying in CPU)', ME.message);
        obj.hCfg.useGPU = 0;
        siteThresh = int16(jrclust.utils.estimateRMS(jrclust.utils.tryGather(samplesIn), 1e5)*obj.hCfg.qqFactor);
    end
end
