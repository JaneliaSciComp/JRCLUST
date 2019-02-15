classdef DetectController < handle
    %DETECTCONTROLLER

    properties (SetAccess=private, SetObservable, Transient)
        hCfg;
        errMsg;
        isError;
    end

    properties (Access=private)
        importTimes; % spike times already known
        importSites; % spike center sites already known

        siteThresh;
        spikeTimes;
        spikeAmps;
        spikeSites;
        centerSites;

        spikesRaw;
        spikesFilt;
        spikeFeatures;
    end

    %% LIFECYCLE
    methods
        function obj = DetectController(hCfg, importTimes, importSites)
            obj.hCfg = hCfg;
            if nargin < 2
                importTimes = [];
            end
            if nargin < 3
                importSites = [];
            end

            obj.importTimes = importTimes(:);
            obj.importSites = importSites(:);

            obj.isError = 0;
        end
    end

    %% USER METHODS
    methods
        function res = detect(obj)
            res = struct();
            t0 = tic();

            % get manually-set spike thresholds
            if ~isempty(obj.hCfg.threshFile)
                try
                    S = load(obj.hCfg.threshFile);
                    siteThresh_ = S.siteThresh;
                    if obj.hCfg.verbose
                        fprintf('Loaded %s\n', obj.hCfg.threshFile);
                    end
                catch ME
                    warning('Could not load threshFile %s: %s', obj.hCfg.threshFile, ME.message);
                    siteThresh_ = [];
                end
            else
                siteThresh_ = [];
            end

            nRecs = numel(obj.hCfg.rawRecordings);
            hRecs = cell(nRecs, 1);

            obj.siteThresh = cell(nRecs, 1);
            obj.spikeTimes = cell(nRecs, 1);
            obj.spikeAmps = cell(nRecs, 1);
            obj.spikeSites = cell(nRecs, 1);
            obj.centerSites = cell(nRecs, 1);
            obj.spikesRaw = cell(nRecs, 1);
            obj.spikesFilt = cell(nRecs, 1);
            obj.spikeFeatures = cell(nRecs, 1);

            recOffset = 0; % sample offset for each recording in sequence

            % load from files
            for iRec = 1:nRecs
                t1 = tic;

                fn = obj.hCfg.rawRecordings{iRec};
                hRec = jrclust.models.recording.Recording(fn, obj.hCfg);

                if hRec.isError
                    error(hRec.errMsg);
                end

                % subset imported samples in this recording interval
                [impTimes, impSites] = deal([]);
                if ~isempty(obj.importTimes)
                    inInterval = (obj.importTimes > recOffset & obj.importTimes <= recOffset + hRec.nSamples);
                    impTimes = obj.importTimes(inInterval) - sampOffset + 1; % shift spike timing

                    % take sites assoc with times between limits
                    if ~isempty(obj.importSites)
                        impSites = obj.importSites(inInterval);
                    end
                end

                recData = detectInRecording(hRec, impTimes, impSites, siteThresh_, obj.hCfg);
                try
                    hRec.setDetections(recData);
                catch ME % maybe rethrow
                    warning('error caught: %s', ME.message);
                    continue;
                end

                t1 = toc(t1);
                nBytesFile = recData.nBytesFile;
                tr = (nBytesFile/jrclust.utils.typeBytes(obj.hCfg.dataType)/obj.hCfg.nChans)/obj.hCfg.sampleRate;

                if obj.hCfg.verbose
                    fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
                       iRec, nRecs, t1, nBytesFile/1e6, nBytesFile/t1/1e6, tr/t1);
                end

                obj.siteThresh{iRec} = recData.siteThresh;
                obj.spikeTimes{iRec} = recData.spikeTimes + recOffset;
                obj.spikeAmps{iRec} = recData.spikeAmps;
                obj.spikeSites{iRec} = recData.spikeSites;
                obj.centerSites{iRec} = recData.centerSites;
                obj.spikesRaw{iRec} = recData.spikesRaw;
                obj.spikesFilt{iRec} = recData.spikesFilt;
                obj.spikeFeatures{iRec} = recData.spikeFeatures;

                recOffset = recOffset + hRec.nSamples;
                hRecs{iRec} = hRec;
            end % for

            res.spikeTimes = cat(1, obj.spikeTimes{:});
            res.spikeAmps = cat(1, obj.spikeAmps{:});
            res.siteThresh = mean(single(cat(1, obj.siteThresh{:})), 1);

            % spike sites
            obj.centerSites = cat(1, obj.centerSites{:});
            res.spikeSites = obj.centerSites(:, 1);
            if size(obj.centerSites, 2) > 1
                res.spikeSites2 = obj.centerSites(:, 2);
            else
                res.spikeSites2 = [];
            end

            % spikes by site
            nSites = obj.hCfg.nSites;
            res.spikesBySite = arrayfun(@(iSite) find(obj.centerSites(:, 1) == iSite), 1:nSites, 'UniformOutput', 0);
            if size(obj.centerSites, 2) >= 2
                res.spikesBySite2 = arrayfun(@(iSite) find(obj.centerSites(:, 2) == iSite), 1:nSites, 'UniformOutput', 0);
            else
                res.spikesBySite2 = cell(1, nSites);
            end
            if size(obj.centerSites, 2) == 3
                res.spikesBySite3 = arrayfun(@(iSite) find(obj.centerSites(:, 3) == iSite), 1:nSites, 'UniformOutput', 0);
            else
                res.spikesBySite3 = [];
            end

            % detected spikes (raw and filtered), features
            res.spikesRaw = cat(3, obj.spikesRaw{:});
            res.rawShape = size(res.spikesRaw);

            res.spikesFilt = cat(3, obj.spikesFilt{:});
            res.filtShape = size(res.spikesFilt);

            res.spikeFeatures = cat(3, obj.spikeFeatures{:});
            res.featuresShape = size(res.spikeFeatures);

            % spike positions
            res.spikePositions = spikePos(res.spikeSites, res.spikeFeatures, obj.hCfg);

            % recordings for inspection
            res.hRecs = hRecs;

            % summarize
            res.detectTime = toc(t0);
            res.detectedOn = now();
        end
    end
end
