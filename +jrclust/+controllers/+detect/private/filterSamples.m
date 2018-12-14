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

        if hCfg.fftThreshMAD > 0
            samplesIn = jrclust.filters.fftClean(samplesIn, hCfg.fftThreshMAD, hCfg.ramToGPUFactor);
        end
    catch ME % GPU denoising failed, retry in CPU
        hCfg.useGPU = false;

        samplesIn = samplesIn_;
        if hCfg.fftThreshMAD > 0
            samplesIn = jrclust.filters.fftClean(samplesIn, hCfg.fftThreshMAD, hCfg.ramToGPUFactor);
        end
    end

    % filter spikes; samples go in padded and come out padded
    try
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], false, hCfg);
    catch ME % GPU filtering failed, retry in CPU
        hCfg.useGPU = false;

        samplesIn = samplesIn_;
        [samplesOut, channelMeans] = jrclust.filters.filtCAR(samplesIn, [], [], false, hCfg);
    end

    clear samplesIn_;

    % common mode rejection
    if hCfg.blank_thresh > 0
        if isempty(channelMeans) % carMode=='whiten'
            channelMeans = jrclust.utils.getCAR(samplesOut, hCfg.carMode, hCfg.ignoreSites);
        end

        keepMe = jrclust.utils.carReject(channelMeans(:), hCfg.blank_period_ms, hCfg.blank_thresh, hCfg.sampleRate);
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