%--------------------------------------------------------------------------
% function [samplesIn, vnWav2_mean] = filtCAR(samplesIn, hCfg, windowPre, windowPost, fTrim_pad)
function [samplesIn, channelMeans] = filtCAR(samplesIn, windowPre, windowPost, trimPad, hCfg)
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

    samplesIn = [windowPre; samplesIn; windowPost];

    % apply filter
    if strcmp(hCfg.filterType, 'user')
        samplesIn = jrclust.filters.userFilter(samplesIn, hCfg.userFiltKernel);
    elseif strcmp(hCfg.filterType, 'fir1')
        samplesIn = jrclust.filters.fir1Filter(samplesIn, ceil(5*hCfg.sampleRate/1000), 2*hCfg.freqLimBP/hCfg.sampleRate);
    elseif strcmp(hCfg.filterType, 'ndiff')
        samplesIn = jrclust.filters.ndiffFilter(samplesIn, hCfg.nDiffOrder);
    elseif strcmp(hCfg.filterType, 'sgdiff')
        samplesIn = jrclust.filters.sgFilter(samplesIn, hCfg.nDiffOrder);
    elseif strcmp(hCfg.filterType, 'bandpass')
        filtOpts = struct('filtOrder', hCfg.filtOrder, ...
                          'sampleRate', hCfg.sampleRate, ...
                          'useElliptic', hCfg.useElliptic, ...
                          'freqLimStop', hCfg.freqLimStop, ...
                          'freqLimBP', hCfg.freqLimBP, ...
                          'freqLimNotch', hCfg.freqLimNotch, ...
                          'nSamplesPad', hCfg.nSamplesPad, ...
                          'useGPUFilt', hCfg.useGPU, ...
                          'gainBoost', hCfg.gainBoost, ...
                          'nDiffOrder', hCfg.nDiffOrder);
        samplesIn = jrclust.filters.bandpassFilter(samplesIn, filtOpts);
    end

    % trim padding
    if trimPad && (nPadPre > 0 || nPadPost > 0)
        samplesIn = samplesIn(nPadPre+1:end-nPadPost, :);
    end

    % global subtraction before
    [samplesIn, channelMeans] = applyCAR(samplesIn, hCfg);
end

%% LOCAL FUNCTIONS
function [samplesIn, channelMeans] = applyCAR(samplesIn, hCfg)
    %APPLYCAR Perform common average referencing (CAR) on filtered traces
    channelMeans = [];

    if strcmp(hCfg.CARMode, 'mean')
        channelMeans = meanExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMeans);
    elseif strcmp(hCfg.CARMode, 'median')
        channelMedians = medianExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMedians);
    end

    samplesIn(:, hCfg.ignoreSites) = 0; % TW do not repair with fMeanSite_drift
end

function means = meanExcluding(samplesIn, ignoreSites)
    %MEANEXCLUDING Calculate mean after excluding ignoreSites
    if isempty(ignoreSites)
        means = int16(mean(samplesIn, 2));
    else
        nSites = size(samplesIn, 2) - numel(ignoreSites);
        means = int16((sum(samplesIn, 2) - sum(samplesIn(:, ignoreSites), 2)) / nSites); % TW BUGFIX
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
