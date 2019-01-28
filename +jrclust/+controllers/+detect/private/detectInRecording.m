function recData = detectInRecording(hRec, impTimes, impSites, presetThresh, hCfg)
    %DETECTINRECORDING Detect spikes in a single Recording
    % get the number of bytes to load from recording
    nBytesFile = subsetBytes(hRec, hCfg.loadTimeLimits);

    % divide recording into several loads, `samples` are columns in data matrix
    [nLoads, nSamplesLoad, nSamplesFinal] = planLoad(nBytesFile, hCfg);

    loadOffset = 0;

    windowPre = [];
    hRec.open();

    siteThresh = cell(nLoads, 1);
    spikeTimes = cell(nLoads, 1);
    spikeAmps = cell(nLoads, 1);
    spikeSites = cell(nLoads, 1);
    centerSites = cell(nLoads, 1);
    spikesRaw = cell(nLoads, 1);
    spikesFilt = cell(nLoads, 1);
    spikeFeatures = cell(nLoads, 1);

    for iLoad = 1:nLoads
        if hCfg.verbose
            t1 = tic;
            fprintf('Processing load %d/%d...\n', iLoad, nLoads);
        end

        if iLoad == nLoads
            nSamples = nSamplesFinal;
        else
            nSamples = nSamplesLoad;
        end

        if hCfg.verbose
            fprintf('\tLoading from file...');
        end

        % load raw samples
        samplesRaw = hRec.readROI(hCfg.siteMap, 1+loadOffset:loadOffset+nSamples);

        % convert samples to int16
        samplesRaw = regularize(samplesRaw, hCfg);
        if hCfg.verbose
            fprintf('done (%0.2f s)\n', toc(t1));
        end

        % load next N samples to ensure we don't miss any spikes at the boundary
        if iLoad < nLoads && hCfg.nSamplesPad > 0
            leftBound = loadOffset + nSamples + 1;
            rightBound = leftBound + hCfg.nSamplesPad - 1;
            windowPost = hRec.readROI(hCfg.siteMap, leftBound:rightBound);
            windowPost = regularize(windowPost, hCfg);
        else
            windowPost = [];
        end

        % denoise and filter samples
        [samplesFilt, keepLower] = filterSamples(samplesRaw, windowPre, windowPost, hCfg);

        % detect spikes in filtered samples
        [intTimes, intSites] = deal([]);
        if ~isempty(impTimes)
            inInterval = (impTimes > loadOffset & impTimes <= loadOffset + nSamples);
            intTimes = impTimes(inInterval) - loadOffset; % shift spike timing

            % take sites assoc with times between limits
            if ~isempty(impSites)
                intSites = impSites(inInterval);
            end
        end

        % detect spikes
        loadData = struct('samplesRaw', samplesRaw, ...
                          'samplesFilt', samplesFilt, ...
                          'keepMe', keepLower, ...
                          'spikeTimes', intTimes, ...
                          'spikeSites', intSites, ...
                          'siteThresh', presetThresh, ...
                          'nPadPre', size(windowPre, 1), ...
                          'nPadPost', size(windowPost, 1));

        % adds spikeTimes, spikeSites; updates siteThresh
        loadData = findPeaks(loadData, hCfg);

        % extract features
        % adds centerSites, spikesRaw, spikesFilt, spikeFeatures
        loadData = extractFeatures(loadData, hCfg);

        % unpack loadData
        siteThresh{iLoad} = loadData.siteThresh;
        spikeTimes{iLoad} = loadData.spikeTimes + loadOffset;
        spikeAmps{iLoad} = loadData.spikeAmps;
        spikeSites{iLoad} = loadData.spikeSites;
        centerSites{iLoad} = loadData.centerSites;
        spikesRaw{iLoad} = loadData.spikesRaw;
        spikesFilt{iLoad} = loadData.spikesFilt;
        spikeFeatures{iLoad} = loadData.spikeFeatures;

        if iLoad < nLoads
            windowPre = samplesRaw(end-hCfg.nSamplesPad+1:end, :);
        end

        clear samplesRaw channelMeans;

        % increment sample offset
        loadOffset = loadOffset + nSamples;
    end % for

    hRec.close();

    recData = struct('siteThresh', cat(1, siteThresh{:}), ...
                     'spikeTimes', cat(1, spikeTimes{:}), ...
                     'spikeAmps', cat(1, spikeAmps{:}), ...
                     'spikeSites', cat(1, spikeSites{:}), ...
                     'centerSites', cat(1, centerSites{:}), ...
                     'spikesRaw', cat(3, spikesRaw{:}), ...
                     'spikesFilt', cat(3, spikesFilt{:}), ...
                     'spikeFeatures', cat(3, spikeFeatures{:}), ...
                     'nBytesFile', nBytesFile);
end

%% LOCAL FUNCTIONS
function nBytesLoad = subsetBytes(hRec, loadTimeLimits)
    %SUBSETBYTES Get number of bytes to load from file, given time limits
    nBytesLoad = hRec.fSizeBytes - hRec.headerOffset;

    if isempty(loadTimeLimits)
        return;
    end

    loadLimits = min(max(loadTimeLimits, 1), hRec.nSamples);
    nSamplesLoad = diff(loadLimits) + 1;
    nBytesLoad = nSamplesLoad * jrclust.utils.typeBytes(hRec.dataType) * hRec.nChans;
end

function samplesRaw = regularize(samplesRaw, hCfg)
    %REGULARIZE Convert samplesRaw to int16
    switch(hCfg.dataType)
        case 'uint16'
            samplesRaw = int16(single(samplesRaw) - 2^15);

        case {'single', 'double'}
            samplesRaw = int16(samplesRaw / hCfg.bitScaling);
    end

    % flip the polarity
    if hCfg.getOr('fInverse_file', 0)
        samplesRaw = -samplesRaw;
    end

    % extract channel means
    if hCfg.tallSkinny
        samplesRaw = samplesRaw';
    else % Catalin's format (samples x channels)
        if ~isempty(hCfg.loadTimeLimits)
            nSamples = size(samplesRaw, 1);
            lims = min(max(round(hCfg.loadTimeLimits * hCfg.sampleRate), 1), nSamples);

            samplesRaw = samplesRaw(lims(1):lims(end), :);
        end
    end
end
