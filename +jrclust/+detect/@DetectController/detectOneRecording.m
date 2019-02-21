function recData = detectOneRecording(obj, hRec, impTimes, impSites, siteThresh)
    %DETECTONERECORDING Detect spikes in a single Recording
    if nargin < 3
        impTimes = [];
    end
    if nargin < 4
        impSites = [];
    end
    if nargin < 5
        siteThresh = [];
    end

    recData = struct('siteThresh', siteThresh(:), ...
                     'spikeTimes', [], ...
                     'spikeAmps', [], ...
                     'spikeSites', [], ...
                     'centerSites', [], ...
                     'spikesRaw', [], ...
                     'spikesFilt', [], ...
                     'spikesFilt2', [], ...
                     'spikesFilt3', [], ...
                     'spikeFeatures', []);

    % divide recording into many loads, samples are columns in data matrix
    [nLoads, nSamplesLoad, nSamplesFinal] = obj.planLoad(hRec);
    loadOffset = 0;
    samplesPre = [];
    hRec.open();
    for iLoad = 1:nLoads
        obj.hCfg.updateLog('procLoad', sprintf('Processing load %d/%d', iLoad, nLoads), 1);

        if iLoad == nLoads
            nSamples = nSamplesFinal;
        else
            nSamples = nSamplesLoad;
        end

        % load raw samples
        iSamplesRaw = hRec.readROI(obj.hCfg.siteMap, 1+loadOffset:loadOffset+nSamples);

        % convert samples to int16
        iSamplesRaw = samplesToInt16(iSamplesRaw, obj.hCfg);

        % load next N samples to ensure we don't miss any spikes at the boundary
        if iLoad < nLoads && obj.hCfg.nSamplesPad > 0
            leftBound = loadOffset + nSamples + 1;
            rightBound = leftBound + obj.hCfg.nSamplesPad - 1;
            samplesPost = hRec.readROI(obj.hCfg.siteMap, leftBound:rightBound);
            samplesPost = samplesToInt16(samplesPost, obj.hCfg);
        else
            samplesPost = [];
        end

        % denoise and filter samples
        obj.hCfg.updateLog('filtSamples', 'Filtering spikes', 1, 0);
        [iSamplesFilt, keepMe] = obj.filterSamples(iSamplesRaw, samplesPre, samplesPost);
        obj.hCfg.updateLog('filtSamples', 'Finished filtering spikes', 0, 1);

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
        loadData = struct('samplesRaw', iSamplesRaw, ...
                          'samplesFilt', iSamplesFilt, ...
                          'keepMe', keepMe, ...
                          'spikeTimes', iSpikeTimes, ...
                          'spikeSites', iSpikeSites, ...
                          'siteThresh', siteThresh, ...
                          'nPadPre', size(samplesPre, 1), ...
                          'nPadPost', size(samplesPost, 1));

        % find peaks: adds spikeAmps, updates spikeTimes, spikeSites,
        %             siteThresh
        loadData = obj.findPeaks(loadData);
        recData.spikeAmps = cat(1, recData.spikeAmps, ...
                                jrclust.utils.tryGather(loadData.spikeAmps));
        recData.spikeSites = cat(1, recData.spikeSites, ...
                                 jrclust.utils.tryGather(loadData.spikeSites));
        recData.siteThresh = cat(1, recData.siteThresh, ...
                                 jrclust.utils.tryGather(loadData.siteThresh));

        % extract spike windows: adds spikesRaw, spikesFilt, centerSites,
        %                        spikesFilt2, spikesFilt3;
        %                        updates spikeTimes
        loadData = obj.samplesToWindows(loadData);
        recData.centerSites = cat(1, recData.centerSites, ...
                                  jrclust.utils.tryGather(loadData.centerSites));
        recData.spikeTimes = cat(1, recData.spikeTimes, ...
                                 jrclust.utils.tryGather(loadData.spikeTimes + loadOffset - size(samplesPre, 1)));
        recData.spikesRaw = cat(3, recData.spikesRaw, ...
                               jrclust.utils.tryGather(loadData.spikesRaw));
        recData.spikesFilt = cat(3, recData.spikesFilt, ...
                                 jrclust.utils.tryGather(loadData.spikesFilt));
        recData.spikesFilt2 = cat(3, recData.spikesFilt2, ...
                                  jrclust.utils.tryGather(loadData.spikesFilt2));
        recData.spikesFilt3 = cat(3, recData.spikesFilt3, ...
                                  jrclust.utils.tryGather(loadData.spikesFilt3));

        % compute features: adds spikeFeatures
        if ~obj.hCfg.getOr('extractAfterDetect', 0)
            loadData = obj.extractFeatures(loadData);
            recData.spikeFeatures = cat(3, recData.spikeFeatures, ...
                                        jrclust.utils.tryGather(loadData.spikeFeatures));
        end

        if iLoad < nLoads
            samplesPre = iSamplesRaw(end-obj.hCfg.nSamplesPad+1:end, :);
        end

        % increment sample offset
        loadOffset = loadOffset + nSamples;

        obj.hCfg.updateLog('procLoad', sprintf('Finished load %d/%d', iLoad, nLoads), 0, 1);
    end % for

    hRec.close();

    if obj.hCfg.getOr('extractAfterDetect', 0) && ~strcmp(obj.hCfg.clusterFeature, 'gpca')
        recData = obj.extractFeatures(recData);
        recData.spikeFeatures = jrclust.utils.tryGather(recData.spikeFeatures);
    end
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
