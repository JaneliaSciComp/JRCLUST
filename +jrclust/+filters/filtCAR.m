function [samplesOut, channelMeans] = filtCAR(samplesIn, windowPre, windowPost, trimPad, hCfg)
    %FILTCAR Apply user-specified filter and common-average referencing
    if nargin < 3
        windowPre = [];
    end
    if nargin < 4
        windowPost = [];
    end
    if nargin < 5
        trimPad = 1;
    end

    nPadPre = size(windowPre, 1);
    nPadPost = size(windowPost, 1);

    samplesOut = [windowPre; samplesIn; windowPost];
    if hCfg.useGPU
        samplesOut = jrclust.utils.tryGpuArray(samplesOut, hCfg.useGPU);
    end

    % apply filter
    if strcmp(hCfg.filterType, 'user')
        samplesOut = jrclust.filters.userFilter(samplesOut, hCfg.userFiltKernel);
    elseif strcmp(hCfg.filterType, 'fir1')
        samplesOut = jrclust.filters.fir1Filter(samplesOut, ceil(5*hCfg.sampleRate/1000), 2*hCfg.freqLimBP/hCfg.sampleRate);
    elseif strcmp(hCfg.filterType, 'ndiff')
        samplesOut = jrclust.filters.ndiffFilter(samplesOut, hCfg.nDiffOrder);
    elseif strcmp(hCfg.filterType, 'sgdiff')
        samplesOut = jrclust.filters.sgFilter(samplesOut, hCfg.nDiffOrder);
    elseif strcmp(hCfg.filterType, 'bandpass')
        hCfg.useGPUFilt = hCfg.useGPU;
        samplesOut = jrclust.filters.bandpassFilter(samplesOut, hCfg);
    end

    % trim padding
    if trimPad && (nPadPre > 0 || nPadPost > 0)
        samplesOut = samplesOut(nPadPre+1:end-nPadPost, :);
    end

    % global subtraction before
    [samplesOut, channelMeans] = applyCAR(samplesOut, hCfg);
    [samplesOut, channelMeans] = jrclust.utils.tryGather(samplesOut, channelMeans);
end

%% LOCAL FUNCTIONS
function [samplesIn, channelMeans] = applyCAR(samplesIn, hCfg)
    %APPLYCAR Perform common average referencing (CAR) on filtered traces
    channelMeans = [];

    if strcmp(hCfg.CARMode, 'mean')
        channelMeans = meanExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMeans);
    elseif strcmp(hCfg.CARMode, 'median')
        channelMeans = medianExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMeans);
    end

    samplesIn(:, hCfg.ignoreSites) = 0; % TW do not repair with fMeanSite_drift
end

function means = meanExcluding(samplesIn, ignoreSites)
    %MEANEXCLUDING Calculate mean after excluding ignoreSites
    if isempty(ignoreSites)
        means = cast(mean(samplesIn, 2), 'like', samplesIn);
    else
        nSites = size(samplesIn, 2) - numel(ignoreSites);
        means = cast((sum(samplesIn, 2) - sum(samplesIn(:, ignoreSites), 2)) / nSites, 'like', samplesIn); % TW BUGFIX
    end
end

function medians = medianExcluding(samplesIn, ignoreSites)
    %MEDIANEXCLUDING Calculate median after excluding ignoreSites
    useGPU = isa(samplesIn, 'gpuArray');
    if useGPU
        samplesIn = jrclust.utils.tryGather(samplesIn);
    end

    if isempty(ignoreSites)
        medians = median(samplesIn, 2);
    else
        goodSites = setdiff(1:size(samplesIn, 2), ignoreSites);
        medians = median(samplesIn(:, goodSites), 2);
    end

    medians = jrclust.utils.tryGpuArray(medians, useGPU);
end
