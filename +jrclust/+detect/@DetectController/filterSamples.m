function [samplesOut, keepMe] = filterSamples(obj, samplesIn, samplesPre, samplesPost)
    %FILTERSAMPLES Denoise and filter raw samples, apply CAR
    samplesIn = [samplesPre; samplesIn; samplesPost];

    if obj.hCfg.fftThresh > 0
        samplesIn = jrclust.filters.fftClean(samplesIn, obj.hCfg.fftThresh, obj.hCfg);
    end

    % filter spikes; samples go in padded and come out padded
    try
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], 0, obj.hCfg);
    catch ME % GPU filtering failed, retry in CPU
        obj.hCfg.updateLog('filtSamples', sprintf('GPU filtering failed: %s (retrying in CPU)', ME.message), 1, 0);

        obj.hCfg.useGPU = 0;
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], 0, obj.hCfg);
    end

    % common mode rejection
    if obj.hCfg.blankThresh > 0
        if isempty(channelMeans) % CARMode=='whiten'
            channelMeans = jrclust.utils.getCAR(samplesOut, obj.hCfg.CARMode, obj.hCfg.ignoreSites);
        end

        keepMe = jrclust.utils.carReject(channelMeans(:), obj.hCfg.blankPeriod, obj.hCfg.blankThresh, obj.hCfg.sampleRate);
        obj.hCfg.updateLog('rejectMotion', sprintf('Rejecting %0.3f %% of time due to motion', (1 - mean(keepMe))*100), 0, 0);
    else
        keepMe = true(size(channelMeans));
    end
end