function recData = detectOneRecording(obj, hRec, fids, impTimes, impSites, siteThresh)
    %DETECTONERECORDING Detect spikes in a single Recording
    if nargin < 4
        impTimes = [];
    end
    if nargin < 5
        impSites = [];
    end
    if nargin < 6
        siteThresh = [];
    end

    rawFid = fids(1);
    filtFid = fids(2);

    recData = struct('siteThresh', siteThresh(:), ...
                     'spikeTimes', [], ...
                     'spikeAmps', [], ...
                     'spikeSites', [], ...
                     'centerSites', [], ...
                     'rawShape', [], ...
                     'filtShape', [], ...
                     'spikesFilt2', [], ...
                     'spikesFilt3', [], ...
                     'spikeFeatures', []);

    % divide recording into many loads, samples are columns in data matrix
    [nLoads, nSamplesLoad, nSamplesFinal] = obj.planLoad(hRec);
    hRec.openRaw();

    % precompute threshold for entire recording (requires 2 passes through data)
    if obj.hCfg.getOr('precomputeThresh', 0)
        samplesPre = [];
        loadOffset = 0;

        obj.hCfg.updateLog('precomputeThresh', 'Precomputing threshold', 1, 0);
        for iLoad = 1:nLoads
            obj.hCfg.updateLog('procLoad', sprintf('Processing load %d/%d', iLoad, nLoads), 1, 0);
            if iLoad == nLoads
                nSamples = nSamplesFinal;
            else
                nSamples = nSamplesLoad;
            end

            % load raw samples
            iSamplesRaw = hRec.readRawROI(obj.hCfg.siteMap, 1+loadOffset:loadOffset+nSamples);

            % convert samples to int16
            iSamplesRaw = samplesToInt16(iSamplesRaw, obj.hCfg);

            % load next N samples to ensure we don't miss any spikes at the boundary
            if iLoad < nLoads && obj.hCfg.nSamplesPad > 0
                leftBound = loadOffset + nSamples + 1;
                rightBound = leftBound + obj.hCfg.nSamplesPad - 1;
                samplesPost = hRec.readRawROI(obj.hCfg.siteMap, leftBound:rightBound);
                samplesPost = samplesToInt16(samplesPost, obj.hCfg);
            else
                samplesPost = [];
            end

            % denoise and filter samples
            obj.hCfg.updateLog('filtSamples', 'Filtering spikes', 1, 0);
            iSamplesRaw = [samplesPre; iSamplesRaw; samplesPost];

            if obj.hCfg.fftThresh > 0
                iSamplesRaw = jrclust.filters.fftClean(iSamplesRaw, obj.hCfg.fftThresh, obj.hCfg);
            end

            % filter spikes; samples go in padded and come out padded
            try
                iSamplesFilt  = jrclust.filters.filtCAR(iSamplesRaw, [], [], 0, obj.hCfg);
            catch ME % GPU filtering failed, retry in CPU
                obj.hCfg.updateLog('filtSamples', sprintf('GPU filtering failed: %s (retrying in CPU)', ME.message), 1, 0);

                obj.hCfg.useGPU = 0;
                iSamplesFilt = jrclust.filters.filtCAR(iSamplesRaw, [], [], 0, obj.hCfg);
            end
%             iSamplesFilt = obj.filterSamples(iSamplesRaw, samplesPre, samplesPost);
            obj.hCfg.updateLog('filtSamples', 'Finished filtering spikes', 0, 1);

            siteThresh = cat(1, siteThresh, obj.computeThreshold(iSamplesFilt));

            nPadPre = size(samplesPre, 1);
            nPadPost = size(samplesPost, 1);
            bounds = [nPadPre + 1, size(iSamplesFilt, 1) - nPadPost]; % inclusive
            success = hRec.writeFilt(iSamplesFilt(bounds(1):bounds(2), :));
            if ~success % abort
                siteThresh = [];
                obj.hCfg.precomputeThresh = 0;
                break;
            end

            if iLoad < nLoads
                samplesPre = iSamplesRaw(end-2*obj.hCfg.nSamplesPad+1:end-obj.hCfg.nSamplesPad, :);
            end

            % increment sample offset
            loadOffset = loadOffset + nSamples;
        end

        siteThresh = mean(siteThresh, 1);

        if isempty(siteThresh)
            obj.hCfg.updateLog('precomputeThresh', 'Failed to precompute threshold', 0, 1);
        else % we made it, don't recompute filtered samples but read them in from disk
            obj.hCfg.updateLog('precomputeThresh', 'Finished precomputing threshold', 0, 1);
            hRec.closeFiltForWriting();
            hRec.openFilt();
        end
    end

    samplesPre = [];
    loadOffset = 0;
    for iLoad = 1:nLoads
        obj.hCfg.updateLog('procLoad', sprintf('Processing load %d/%d', iLoad, nLoads), 1);

        if iLoad == nLoads
            nSamples = nSamplesFinal;
        else
            nSamples = nSamplesLoad;
        end

        % load raw samples
        iSamplesRaw = hRec.readRawROI(obj.hCfg.siteMap, 1+loadOffset:loadOffset+nSamples);

        % convert samples to int16
        iSamplesRaw = samplesToInt16(iSamplesRaw, obj.hCfg);

        % load next N samples to ensure we don't miss any spikes at the boundary
        if iLoad < nLoads && obj.hCfg.nSamplesPad > 0
            leftBound = loadOffset + nSamples + 1;
            rightBound = leftBound + obj.hCfg.nSamplesPad - 1;
            samplesPost = hRec.readRawROI(obj.hCfg.siteMap, leftBound:rightBound);
            samplesPost = samplesToInt16(samplesPost, obj.hCfg);
        else
            samplesPost = [];
        end

        nPadPre = size(samplesPre, 1);
        nPadPost = size(samplesPost, 1);

        % samples pre-filtered, read from disk
        if obj.hCfg.getOr('precomputeThresh', 0)
            obj.hCfg.updateLog('filtSamples', 'Reading filtered samples from disk', 1, 0);
            iSamplesFilt = hRec.readFiltROI(1:obj.hCfg.nSites, 1+loadOffset-nPadPre:loadOffset+nSamples+nPadPost)';
            obj.hCfg.updateLog('filtSamples', sprintf('Read %d filtered samples', size(iSamplesFilt, 1)), 0, 1);

            % common mode rejection
            if obj.hCfg.blankThresh > 0
                channelMeans = jrclust.utils.getCAR(iSamplesFilt, obj.hCfg.CARMode, obj.hCfg.ignoreSites);

                keepMe = jrclust.utils.carReject(channelMeans(:), obj.hCfg.blankPeriod, obj.hCfg.blankThresh, obj.hCfg.sampleRate);
                obj.hCfg.updateLog('rejectMotion', sprintf('Rejecting %0.3f %% of time due to motion', (1 - mean(keepMe))*100), 0, 0);
            else
                keepMe = true(size(iSamplesFilt, 1), 1);
            end
        else
            % denoise and filter samples
            obj.hCfg.updateLog('filtSamples', 'Filtering samples', 1, 0);
            [iSamplesFilt, keepMe] = obj.filterSamples(iSamplesRaw, samplesPre, samplesPost);
            obj.hCfg.updateLog('filtSamples', 'Finished filtering samples', 0, 1);
        end

        % detect spikes in filtered samples
        [iSpikeTimes, iSpikeSites] = deal([]);
        if ~isempty(impTimes)
            inInterval = (impTimes > loadOffset & impTimes <= loadOffset + nSamples);
            iSpikeTimes = impTimes(inInterval) - loadOffset; % shift spike timing

            % take sites assoc with times between limits
            if ~isempty(impSites)
                iSpikeSites = impSites(inInterval);
            end
        end

        %%% detect spikes
        loadData = struct('samplesRaw', [samplesPre; iSamplesRaw; samplesPost], ...
                          'samplesFilt', iSamplesFilt, ...
                          'keepMe', keepMe, ...
                          'spikeTimes', iSpikeTimes, ...
                          'spikeSites', iSpikeSites, ...
                          'siteThresh', siteThresh, ...
                          'nPadPre', nPadPre, ...
                          'nPadPost', nPadPost);

        % find peaks: adds spikeAmps, updates spikeTimes, spikeSites,
        %             siteThresh
        loadData = obj.findPeaks(loadData);
        recData.spikeAmps = cat(1, recData.spikeAmps, loadData.spikeAmps);
        recData.spikeSites = cat(1, recData.spikeSites, loadData.spikeSites);
        recData.siteThresh = [recData.siteThresh, loadData.siteThresh];

        % extract spike windows: adds centerSites, updates spikeTimes
        loadData = obj.samplesToWindows(loadData);
        recData.centerSites = cat(1, recData.centerSites, loadData.centerSites);
        recData.spikeTimes = cat(1, recData.spikeTimes, loadData.spikeTimes + loadOffset - size(samplesPre, 1));
        recData.spikesFilt2 = cat(3, recData.spikesFilt2, loadData.spikesFilt2);
        recData.spikesFilt3 = cat(3, recData.spikesFilt3, loadData.spikesFilt3);

        % write out spikesRaw and update shape
        fwrite(rawFid, loadData.spikesRaw, '*int16');
        if isempty(recData.rawShape)
            recData.rawShape = size(loadData.spikesRaw);
        else
            recData.rawShape(3) = recData.rawShape(3) + size(loadData.spikesRaw, 3);
        end

        % write out spikesFilt and update shape
        fwrite(filtFid, loadData.spikesFilt, '*int16');
        if isempty(recData.filtShape)
            recData.filtShape = size(loadData.spikesFilt);
        else
            recData.filtShape(3) = recData.filtShape(3) + size(loadData.spikesFilt, 3);
        end

        % compute features: adds spikeFeatures
        if ~obj.hCfg.getOr('extractAfterDetect', 0)
            loadData = obj.extractFeatures(loadData);
            recData.spikeFeatures = cat(3, recData.spikeFeatures, loadData.spikeFeatures);
        end

        if iLoad < nLoads
            samplesPre = iSamplesRaw(end-obj.hCfg.nSamplesPad+1:end, :);
        end

        % increment sample offset
        loadOffset = loadOffset + nSamples;

        jrclust.utils.tryGather(iSamplesFilt);
        obj.hCfg.updateLog('procLoad', sprintf('Finished load %d/%d', iLoad, nLoads), 0, 1);
    end % for

    hRec.closeRaw();
    hRec.closeFilt();
end

%% LOCAL FUNCTIONS
function samplesRaw = samplesToInt16(samplesRaw, hCfg)
    %SAMPLESTOINT16 Convert samplesRaw to int16
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
