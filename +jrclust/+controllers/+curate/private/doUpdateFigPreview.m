%function hFigPreview = doUpdateFigPreview(hFigPreview, S_fig, fKeepView)
function hFigPreview = doUpdateFigPreview(hFigPreview, figData, fKeepView, hCfg)
    %DOUPDATEFIGPREVIEW Handle parameter update requests
    hFigPreview.wait(true);
    fftThreshMAD = figData.fftThreshMAD;

    if fftThreshMAD > 0
        figData.tracesClean = jrclust.filters.fftClean(figData.tracesRaw, figData.fftThreshMAD, hCfg.ramToGPUFactor); % fft filter
    else
        figData.tracesClean = figData.tracesRaw;
    end

    % Find bad sites
    if figData.siteCorrThresh > 0
        figData.ignoreMe = (figData.maxCorrSite < figData.siteCorrThresh);
        figData.ignoreSites = find(figData.ignoreMe);
    elseif ~isempty(figData.ignoreSites)
        figData.ignoreSites = figData.ignoreSites;
        figData.ignoreMe = false(size(figData.maxCorrSite));
        figData.ignoreMe(figData.ignoreSites) = 1;
        figData.ignoreSites = [];
    else
        figData.ignoreMe = false(size(figData.maxCorrSite));
        figData.ignoreSites = [];
    end

    % Perform filter fftThreshMAD
    hCfg.setTemporaryParams('carMode', 'none', 'useGPU', 0, 'filterType', figData.filterType, ...
        'blank_period_ms', figData.blank_period_ms, 'blankThresh', figData.blankThresh, 'fParfor', 0);

    tracesFilt = jrclust.filters.filtCAR(figData.tracesClean, [], [], false, hCfg);
    tracesCAR = jrclust.utils.getCAR(tracesFilt, figData.carMode, figData.ignoreSites);

    if ~strcmpi(figData.carMode, 'none')
        tracesFilt = bsxfun(@minus, tracesFilt, int16(tracesCAR));
    end

    figData.tracesCAR = jrclust.utils.madScore(mean(tracesCAR, 2)); % Save in MAD unit
    figData.tracesFilt = tracesFilt;

    % reference threshold
    [rawPSD, figData.vrFreq_psd] = getPSD(figData.tracesRaw(:, ~figData.ignoreMe), hCfg.sampleRate, 4);
    cleanPSD = getPSD(figData.tracesClean(:, ~figData.ignoreMe), hCfg.sampleRate, 4);

    figData.vrPower_psd = mean(rawPSD, 2);
    figData.vrPower_clean_psd = mean(cleanPSD, 2);

    % Apply threshold and perform spike detection
    siteRMS = jrclust.utils.estimateRMS(tracesFilt, 1e5);
    siteThresh = int16(siteRMS * figData.qqFactor);
    siteThresh(figData.ignoreMe) = 0;
    figData.siteThresh = siteThresh;

    figData.isThreshCrossing = bsxfun(@lt, tracesFilt, -abs(siteThresh));
    figData.isThreshCrossing(:, figData.ignoreMe) = false; % ignore threshold crossings on bad sites

    % Spike detection
    [keepMe, figData.channelMeansMAD] = jrclust.utils.carReject(tracesCAR, hCfg.blank_period_ms, hCfg.blankThresh, hCfg.sampleRate);
    [figData.spikeTimes, figData.spikeAmps, spikeSites] = jrclust.utils.detectPeaks(tracesFilt, siteThresh, keepMe, hCfg);

    durationSecs = size(tracesFilt, 1) / hCfg.sampleRate;

    figData.siteEventRate = hist(spikeSites, 1:hCfg.nSites)/durationSecs; % event count
    figData.siteEventSNR = abs(single(arrayfun(@(i) median(figData.spikeAmps(spikeSites == i)), 1:hCfg.nSites)))./siteRMS;

    figData.keepMe = keepMe;
    figData.spikeSites = spikeSites;
    % Spike stats: such as # sites/event over threshold

    hCfg.resetTemporaryParams();

    hFigPreview.figData = figData;
    hFigPreview = doPlotFigPreview(hFigPreview, hCfg, fKeepView);
    hFigPreview.wait(false);
end

%% LOCAL FUNCTIONS
function [mrPower, vrFreq] = getPSD(traces, sampleRate, nSkip)
    %GETPSD Compute power spectral density of traces
    if nSkip > 1
        traces = traces(1:nSkip:end, :);
    end

    n = size(traces, 1);
    n1 = round(n/2);

    mrPower = fft(jrclust.utils.meanSubtract(single(traces)));
    mrPower = jrclust.utils.pow2db(abs(mrPower(2:n1 + 1, :))) / n;

    if nargout > 1
        vrFreq = sampleRate*(1:n1)'/n;
    end
end

