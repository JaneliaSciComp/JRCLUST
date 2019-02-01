function [samplesOut, keepMe] = filterSamples(samplesIn, windowPre, windowPost, hCfg)
    %FILTERSAMPLES Denoise and filter raw samples, apply CAR
    tFilt = tic;

    if hCfg.verbose
        fprintf('\tFiltering spikes...');
    end
    samplesIn = [windowPre; samplesIn; windowPost];

    samplesIn_ = samplesIn; % keep a copy in CPU
    try
        samplesIn = jrclust.utils.tryGpuArray(samplesIn, hCfg.useGPU);

        if hCfg.fftThresh > 0
            samplesIn = jrclust.filters.fftClean(samplesIn, hCfg.fftThresh, hCfg.ramToGPUFactor);
        end
    catch ME % GPU denoising failed, retry in CPU
        hCfg.useGPU = 0;

        samplesIn = samplesIn_;
        if hCfg.fftThresh > 0
            samplesIn = jrclust.filters.fftClean(samplesIn, hCfg.fftThresh, hCfg.ramToGPUFactor);
        end
    end

    % filter spikes; samples go in padded and come out padded
    try
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], 0, hCfg);
    catch ME % GPU filtering failed, retry in CPU
        hCfg.useGPU = 0;

        samplesIn = samplesIn_;
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], 0, hCfg);
    end

    clear samplesIn_;

    % common mode rejection
    if hCfg.blankThresh > 0
        if isempty(channelMeans) % CARMode=='whiten'
            channelMeans = jrclust.utils.getCAR(samplesOut, hCfg.CARMode, hCfg.ignoreSites);
        end

        keepMe = jrclust.utils.carReject(channelMeans(:), hCfg.blankPeriod, hCfg.blankThresh, hCfg.sampleRate);
        if hCfg.verbose
            fprintf('!! rejecting %0.3f %% of time due to motion !!', (1 - mean(keepMe))*100);
        end
    else
        keepMe = true(size(channelMeans));
    end

    if hCfg.verbose
        fprintf('\tdone (%0.2f) s\n', toc(tFilt));
    end
end