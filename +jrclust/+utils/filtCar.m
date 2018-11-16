%--------------------------------------------------------------------------
% function [samplesIn, vnWav2_mean] = filtCar(samplesIn, hCfg, windowPre, windowPost, fTrim_pad)
function [samplesIn, channelMeans] = filtCar(samplesIn, windowPre, windowPost, trimPad, hCfg)
    %FILTCAR Apply user-specified filter and common-average referencing

    if nargin < 3
        windowPre = [];
    end
    if nargin < 4
        windowPost = [];
    end
    if nargin < 5
        trimPad = true;
    end

    nPadPre = size(windowPre, 1);
    nPadPost = size(windowPost, 1);

    samplesIn = [windowPre; samplesIn; windowPost];

    % apply filter
    if strcmp(hCfg.filterType, 'user')
        samplesIn = jrclust.filters.userFilter(samplesIn, hCfg.userFiltKernel);
    elseif strcmp(hCfg.filterType, 'fir1')
        samplesIn = jrclust.filters.fir1Filter(samplesIn, ceil(5*hCfg.sampleRate/1000), 2*hCfg.freqLim/hCfg.sampleRate);
    elseif strcmp(hCfg.filterType, 'ndiff')
        samplesIn = jrclust.filters.ndiffFilter(samplesIn, hCfg.nDiff_filt);
    elseif strcmp(hCfg.filterType, 'fftdiff')
        samplesIn = jrclust.filters.fftdiffFilter(samplesIn, 2*hCfg.freqLim/hCfg.sampleRate, hCfg.ramToGPUFactor);
    elseif any(strcmp(hCfg.filterType, {'sgdiff', 'sgfilt'}))
        samplesIn = jrclust.filters.sgFilter(samplesIn, hCfg.nDiff_filt);
    elseif strcmp(hCfg.filterType, 'bandpass')
        filtOpts = struct('filtOrder', hCfg.filtOrder, ...
                          'sampleRate', hCfg.sampleRate, ...
                          'useElliptic', hCfg.useElliptic, ...
                          'freqLimStop', hCfg.freqLimStop, ...
                          'freqLim', hCfg.freqLim, ...
                          'freqLimNotch', hCfg.freqLimNotch, ...
                          'nSamplesPad', hCfg.nSamplesPad, ...
                          'useGPUFilt', hCfg.useGPU, ...
                          'gainBoost', hCfg.gainBoost, ...
                          'nDiff_filt', hCfg.nDiff_filt);
        samplesIn = jrclust.filters.bandpassFilter(samplesIn, filtOpts);
    elseif strcmp(hCfg.filterType, 'ndist')
        samplesIn = jrclust.filters.ndistFilter(samplesIn, hCfg.ndist_filt);
    elseif ~any(strcmp(hCfg.filterType, {'none', 'skip'}))
        error('invalid filter option (filterType=''%s'')', filterType);
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

    if any(strcmpi(hCfg.carMode, {'tmean', 'nmean'})) % dead code (deprecated option?)
        trimLim = [.25, .75];
        miSite_ref = jrclust.utils.findSiteNeighbors(hCfg.mrSiteXY, hCfg.nSites_ref + hCfg.nSites_excl_ref, hCfg.ignoreSites, hCfg.viShank_site);
        miSite_ref = miSite_ref(hCfg.nSites_excl_ref+1:end, :); %excl three nearest sites
        viChan_keep = round(trimLim * size(miSite_ref,1));
        viChan_keep = (viChan_keep(1)+1):viChan_keep(2);
        mnWav1_pre = samplesIn;

        if strcmpi(hCfg.carMode, 'tmean')
            for iChan = 1:size(samplesIn,2)
                mnWav2 = sort(mnWav1_pre(:, miSite_ref(:,iChan)), 2);
                gvr_tmean = sum(mnWav2(:, viChan_keep), 2); %may go out of range
                gvr_tmean = int16(single(gvr_tmean)/numel(viChan_keep));
                samplesIn(:, iChan) = mnWav1_pre(:,iChan) - gvr_tmean;
                fprintf('.');
            end
        else
            for iChan = 1:size(samplesIn,2)
                gvr_tmean = sum(mnWav1_pre(:, miSite_ref(:,iChan)), 2); %may go out of range
                gvr_tmean = int16(single(gvr_tmean)/size(miSite_ref,1));
                samplesIn(:, iChan) = mnWav1_pre(:, iChan) - gvr_tmean;
                fprintf('.');
            end
        end
    elseif strcmpi(hCfg.carMode, 'mean')
        channelMeans = meanExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMeans);
    elseif strcmpi(hCfg.carMode, 'median')
        channelMedians = medianExcluding(samplesIn, hCfg.ignoreSites);
        samplesIn = bsxfun(@minus, samplesIn, channelMedians);
    elseif strcmpi(hCfg.carMode, 'whiten')  
        samplesIn = whiten(samplesIn, hCfg.ignoreSites, hCfg.ramToGPUFactor);
    end

    samplesIn(:, hCfg.ignoreSites) = 0; % TW do not repair with fMeanSite_drift

    % TODO: evaluate this code
    % % ignoreSites should be treated carefully. try to repair using nearest sites?
    % if get_(hCfg, 'fMeanSite_drift')
    %     mnWav1 = meanSite_drift_(mnWav1, hCfg);
    % elseif fRepairSites
    %     mnWav1 = meanSite_drift_(mnWav1, hCfg, hCfg.ignoreSites);
    % else
    %     mnWav1(:, hCfg.ignoreSites) = 0;
    % end
    % fprintf('\n\ttook %0.1fs.\n', toc(t1));
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

function samplesOut = whiten(samplesIn, ignoreSites, ramToGPUFactor)
    %WHITEN Apply spatial whitening to samplesIn
    nSamplesMax = round(size(samplesIn, 1)/ramToGPUFactor);

    fprintf('Whitening\n\t');
    tw = tic;
    
    % TODO: address subsample_mr_
    subsamplesIn = subsample_mr_(samplesIn, nSamplesMax, 1);

    goodSites = setdiff(1:size(samplesIn, 2), ignoreSites);
    if ~isempty(ignoreSites)
        subsamplesIn = subsamplesIn(:, goodSites);
    end
    subsamplesIn = single(subsamplesIn);

    mrXXT = subsamplesIn' * subsamplesIn;
    [U,D] = eig(mrXXT + eps('single'));
    Sinv = diag(1./sqrt(diag(D)));

    scale = mean(sqrt(diag(mrXXT)));
    wmat = (U * Sinv * U') * scale; % whitening matrix

    % apply whitening matrix
    samplesOut = zeros(size(samplesIn), 'like', samplesIn);
    if ~isempty(ignoreSites)
        samplesIn = samplesIn(:, goodSites);
    end

    samplesIn = single(samplesIn);
    for iSite_ = 1:numel(goodSites)
        iSite = goodSites(iSite_);
        samplesOut(:, iSite) = int16(samplesIn * wmat(:, iSite_));
        fprintf('.');
    end

    fprintf('\n\ttook %0.1fs\n', toc(tw));
end


