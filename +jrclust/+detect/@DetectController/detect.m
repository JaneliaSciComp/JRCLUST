function res = detect(obj)
    %DETECT Detect spikes in all recordings
    tDetect = tic();

    % get manually-set spike thresholds
    if ~isempty(obj.hCfg.threshFile)
        try
            S = load(obj.hCfg.threshFile);
            siteThresh = S.siteThresh;
            obj.hCfg.updateLog('loadThresh', sprintf('Loaded threshFile %s', obj.hCfg.threshFile), 0, 0);
        catch ME
            warning('Could not load threshFile %s: %s', obj.hCfg.threshFile, ME.message);
            siteThresh = [];
        end
    else
        siteThresh = [];
    end

    nRecs = numel(obj.hCfg.rawRecordings);
    hRecs = cell(nRecs, 1);

    res = struct('siteThresh', siteThresh, ...
                 'spikeTimes', [], ...
                 'spikeAmps', [], ...
                 'centerSites', [], ...
                 'rawShape', [], ...
                 'filtShape', [], ...
                 'spikesFilt2', [], ...
                 'spikesFilt3', [], ...
                 'spikeFeatures', []);

    recOffset = 0; % sample offset for each recording in sequence

    rawFid = fopen(obj.hCfg.rawFile, 'w');
    filtFid = fopen(obj.hCfg.filtFile, 'w');

    % load from files
    for iRec = 1:nRecs
        if nRecs > 1 % reset random seeds for each file for reproducibility
            if obj.hCfg.useGPU
                parallel.gpu.rng(obj.hCfg.randomSeed);
            end
            rng(obj.hCfg.randomSeed);
        end

        fn = obj.hCfg.rawRecordings{iRec};
        hRec = jrclust.detect.newRecording(fn, obj.hCfg);

        if hRec.isError
            error(hRec.errMsg);
        end

        % subset imported samples in this recording interval
        [impTimes, impSites] = deal([]);
        if ~isempty(obj.importTimes)
            inInterval = (obj.importTimes > recOffset & obj.importTimes <= recOffset + hRec.nSamples);
            impTimes = obj.importTimes(inInterval) - recOffset; % shift spike timing

            % take sites assoc with times between limits
            if ~isempty(obj.importSites)
                impSites = obj.importSites(inInterval);
            end
        end

        obj.hCfg.updateLog('fileLoad', sprintf('Processing file %s (%d/%d)', hRec.rawPath, iRec, nRecs), 1, 0);
        recData = obj.detectOneRecording(hRec, [rawFid, filtFid], impTimes, impSites, siteThresh);

        obj.hCfg.updateLog('fileLoad', sprintf('Finished processing file %s (%d/%d)', hRec.rawPath, iRec, nRecs), 0, 1);

        res.siteThresh = [res.siteThresh, recData.siteThresh];
        res.spikeTimes = cat(1, res.spikeTimes, recData.spikeTimes + recOffset);
        res.spikeAmps = cat(1, res.spikeAmps, recData.spikeAmps);
        res.centerSites = cat(1, res.centerSites, recData.centerSites);
        res.spikesFilt2 = cat(3, res.spikesFilt2, recData.spikesFilt2);
        res.spikesFilt3 = cat(3, res.spikesFilt3, recData.spikesFilt3);
        if isfield(recData, 'spikeFeatures')
            res.spikeFeatures = cat(3, res.spikeFeatures, recData.spikeFeatures);
        end

        % update rawShape, filtShape
        if isempty(res.rawShape)
            res.rawShape = recData.rawShape;
        else
            res.rawShape(3) = res.rawShape(3) + size(recData.rawShape, 3);
        end
        if isempty(res.filtShape)
            res.filtShape = recData.filtShape;
        else
            res.filtShape(3) = res.filtShape(3) + size(recData.filtShape, 3);
        end

        recOffset = recOffset + hRec.nSamples;
        hRecs{iRec} = hRec;
    end % for

    % close raw/filtered spike files
    try
        fclose(rawFid);
        fclose(filtFid);
    catch ME
    end

    % load in spikesRaw, spikesFilt
    fid = fopen(obj.hCfg.rawFile, 'r');
    res.spikesRaw = reshape(fread(fid, inf, '*int16'), res.rawShape);
    fclose(fid);
    fid = fopen(obj.hCfg.filtFile, 'r');
    res.spikesFilt = reshape(fread(fid, inf, '*int16'), res.filtShape);
    fclose(fid);

    % compute features from all spikes over all recordings
    if obj.hCfg.getOr('extractAfterDetect', 0) || strcmp(obj.hCfg.clusterFeature, 'gpca')
        res = obj.extractFeatures(res);
    end

    featuresFid = fopen(obj.hCfg.featuresFile, 'w');
    fwrite(featuresFid, res.spikeFeatures, '*single');
    fclose(featuresFid);

    % compute the mean of the siteThresh from each recording
    res.meanSiteThresh = mean(single(res.siteThresh), 2);

    % spike sites
    res.spikeSites = res.centerSites(:, 1);
    if size(res.centerSites, 2) > 1
        res.spikeSites2 = res.centerSites(:, 2);
    else
        res.spikeSites2 = [];
    end
    if size(res.centerSites, 2) > 2
        res.spikeSites3 = res.centerSites(:, 3);
    else
        res.spikeSites3 = [];
    end

    % spikes by site
    nSites = obj.hCfg.nSites;
    res.spikesBySite = arrayfun(@(iSite) find(res.centerSites(:, 1) == iSite), 1:nSites, 'UniformOutput', 0);
    if size(res.centerSites, 2) >= 2
        res.spikesBySite2 = arrayfun(@(iSite) find(res.centerSites(:, 2) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        res.spikesBySite2 = cell(1, nSites);
    end
    if size(res.centerSites, 2) == 3
        res.spikesBySite3 = arrayfun(@(iSite) find(res.centerSites(:, 3) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        res.spikesBySite3 = cell(1, nSites);
    end

    % detected spikes (raw and filtered), features
    res.featuresShape = size(res.spikeFeatures);

    % spike positions
    res.spikePositions = obj.spikePos(res.spikeSites, res.spikeFeatures);

    % recordings for inspection
    res.hRecs = hRecs;

    % summarize
    res.detectTime = toc(tDetect);
    res.detectedOn = now();
end
