function samplesOut = filtTouV(samplesIn, hCfg)
    %FILTTOUV Transform mean-subtracted filtered samples to units of uV
    samplesOut = jrclust.utils.bit2uV(samplesIn, hCfg);
    samplesOut = jrclust.utils.meanSubtract(samplesOut);
end
