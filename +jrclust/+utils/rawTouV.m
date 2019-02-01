function samplesOut = rawTouV(samplesIn, hCfg)
    %RAWTOUV Transform mean-subtracted raw samples to units of uV
    samplesOut = jrclust.utils.meanSubtract(single(samplesIn) * hCfg.bitScaling);
end
