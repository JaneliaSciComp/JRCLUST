classdef DetectionController < handle
    %DETECTIONCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        hCfg;
        hRecs;
        isError;
    end

    properties
        knownTimes; % spike times already known
        knownSites; % spike center sites already known

        siteThresh; % site-wise detection thresholds
        spikeTimes; % detected spike times
        spikeAmps;  % detected spike amplitudes
        spikeSites; % detected spike center sites

        spikesRaw;  % raw waveform tensor, nSamples x nChannels x nSpikes
        spikesFilt; % filtered waveform tensor, nSamples x nChannels x nSpikes
        spikeFeatures; % spike features tensor
        miSite_spk; % 
    end

    % LIFECYCLE
    methods
        function obj = DetectionController(hCfg, knownTimes, knownSites)
            obj.hCfg = hCfg;
            obj.hRecs = jrclust.models.Recording.empty;
            if nargin < 2
                knownTimes = [];
            end
            if nargin < 3
                knownSites = [];
            end

            obj.knownTimes = knownTimes(:);
            obj.knownSites = knownSites(:);

            obj.siteThresh = {};
            obj.spikeTimes = {};
            obj.spikeAmps = {};
            obj.spikeSites = {};
            obj.miSite_spk = {};
            obj.spikesRaw = {};
            obj.spikesFilt = {};
            obj.spikeFeatures = {};

            obj.isError = false;
        end
    end

    % USER METHODS
    methods
        function res = detect(obj)
            res = struct();
            t0 = tic;

            dtype = obj.hCfg.dtype;
            nChans = obj.hCfg.nChans;
            headerOffset = obj.hCfg.headerOffset;

            % get manually-set spike thresholds
            if ~isempty(obj.hCfg.threshFile)
                try
                    S = load(obj.hCfg.threshFile);
                    siteThresh_ = S.siteThresh;
                    fprintf('Loaded %s\n', obj.hCfg.threshFile);
                catch
                    warning('Could not load threshFile %s', obj.hCfg.threshFile);
                    siteThresh_ = [];
                end
            else
                siteThresh_ = [];
            end

            nRecs = numel(obj.hCfg.rawRecordings);

            % load from files
            for iRec = 1:nRecs
                t1 = tic;

                fn = obj.hCfg.rawRecordings{iRec};
                rec = jrclust.models.Recording(fn, dtype, nChans, headerOffset);

                % get the number of bytes to load from recording
                nBytesFile = rec.subsetBytes(obj.hCfg.loadTimeLimits);
                % divide recording into many loads, `samples` are columns in data matrix
                [nLoads, nSamplesLoad, nSamplesFinal] = obj.planLoad(nBytesFile);

                % current column in data matrix
                sampOffset = 1;

                windowPre = [];
                for iLoad = 1:nLoads
                    t1 = tic;

                    fprintf('Processing %d/%d of file %d/%d...\n', iLoad, nLoads, iRec, nRecs);
                    if iLoad == nLoads
                        nSamples = nSamplesFinal;
                    else
                        nSamples = nSamplesLoad;
                    end
                    fprintf('\tLoading from file...');

                    % load raw samples
                    samplesRaw = rec.readROI(obj.hCfg.siteMap, sampOffset:sampOffset+nSamples-1);

                    % convert samples to int16 and get channel means
                    [samplesRaw, channelMeans] = obj.regularize(samplesRaw);
                    fprintf('done (%0.2f s)\n', toc(t1));

                    % select given spikes in this time interval
                    [intTimes, intSites] = deal([]);
                    if ~isempty(obj.knownTimes)
                        inInterval = (obj.knownTimes >= sampOffset & obj.knownTimes < sampOffset + nSamples);
                        intTimes = obj.knownTimes(inInterval) - sampOffset + 1; % shift spike timing

                        % take sites assoc with times between limits
                        if ~isempty(obj.knownSites)
                            intSites = obj.knownSites(inInterval);
                        end
                    end

                    % load next N samples to ensure we don't miss any spikes at the boundary
                    if iLoad < nLoads && obj.hCfg.nSamplesPad > 0
                        windowPost = rec.readROI(obj.hCfg.siteMap, sampOffset:sampOffset+nSamples+obj.hCfg.nSamplesPad-1);
                        windowPost = obj.regularize(windowPost);
                    else
                        windowPost = [];
                    end

                    % denoise and filter samples
                    [samplesFilt, keepMe] = obj.filterSamples(samplesRaw, windowPre, windowPost);
                    % switch obj.hCfg.vcFilter_detect (DEPRECATED?)
                    %     case {'', 'none'}
                    %         mnWav3 = samplesFilt;

                    %     case 'ndist'
                    %         [mnWav3, nShift_post] = filter_detect_(samplesRaw, obj.hCfg); % pass raw trace

                    %     otherwise
                    %         [mnWav3, nShift_post] = filter_detect_(samplesFilt, obj.hCfg); % pass filtered trace
                    % end

                    % detect spikes in filtered samples
                    obj.findSpikes(samplesFilt, keepMe, intTimes, intSites, siteThresh_, size(windowPre, 1), size(windowPost, 1));

                    % extract features
                    obj.extractFeatures(samplesRaw, samplesFilt, size(windowPre, 1));

                    obj.spikeTimes{end} = obj.spikeTimes{end} + sampOffset - 1;

                    % increment sample offset
                    sampOffset = sampOffset + nSamples;

                    if iLoad < nLoads
                        windowPre = samplesRaw(end-obj.hCfg.nSamplesPad+1:end, :);
                    end

                    clear samplesRaw channelMeans;
                end
                t1 = toc(t1);
                tr = (nBytesFile/jrclust.utils.typeBytes(obj.hCfg.dtype)/obj.hCfg.nChans)/obj.hCfg.sampleRate;
                fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
                    iRec, nRecs, t1, nBytesFile/1e6, nBytesFile/t1/1e6, tr/t1);
            end

            res.runtime = toc(t0);
            obj.miSite_spk = cat(1, obj.miSite_spk{:});
            obj.spikeTimes = cat(1, obj.spikeTimes{:});
            obj.spikeAmps = cat(1, obj.spikeAmps{:});
            obj.siteThresh = cat(1, obj.siteThresh{:});
        end
    end

    % UTILITY METHODS
    methods (Access=protected, Hidden)
        function [nLoads, nSamplesLoad, nSamplesFinal] = planLoad(obj, nBytesFile)
            %PLANLOAD Get number of samples to load in each chunk of a file
            bps = jrclust.utils.typeBytes(obj.hCfg.dtype);

            % nColumns in data matrix
            nSamples = floor(nBytesFile / bps / obj.hCfg.nChans);

            % if not constrained by user, try to compute maximum bytes/load
            if isempty(obj.hCfg.maxBytesLoad)
                if obj.hCfg.useGPU
                    S = gpuDevice(); % does not reset GPU
                    nBytes = floor(S(1).AvailableMemory());
                elseif ispc()
                    S = memory();
                    nBytes = floor(S.MaxPossibleArrayBytes());
                else % no hints given, assume infinite memory
                    nBytes = inf;
                end

                obj.hCfg.maxBytesLoad = floor(nBytes / obj.hCfg.gpuLoadFactor);
            end

            % if not constrained by user, try to compute maximum samples/load
            if isempty(obj.hCfg.maxSecLoad)
                nSamplesMax = floor(obj.hCfg.maxBytesLoad / obj.hCfg.nChans / bps);
            else
                nSamplesMax = floor(obj.hCfg.sampleRate * obj.hCfg.maxSecLoad);
            end

            if ~obj.hCfg.fTranspose_bin % load entire file, Catalin's format
                [nLoads, nSamplesLoad, nSamplesFinal] = deal(1, nSamples, nSamples);
            else
                [nLoads, nSamplesLoad, nSamplesFinal] = jrclust.utils.partitionLoad(nSamples, nSamplesMax);
            end
        end

        function [samplesRaw, channelMeans] = regularize(obj, samplesRaw)
            %REGULARIZE Convert samplesRaw to int16 and compute channel means
            switch(obj.hCfg.dtype)
                case 'uint16'
                    samplesRaw = int16(single(samplesRaw) - 2^15);

                case {'single', 'double'}
                    samplesRaw = int16(samplesRaw / obj.hCfg.bitScaling);
            end

            % flip the polarity
            if obj.hCfg.fInverse_file
                samplesRaw = -samplesRaw;
            end

            % extract channel means
            if obj.hCfg.fTranspose_bin
                if nargout > 1
                    channelMeans = single(mean(samplesRaw, 1));
                end
                samplesRaw = samplesRaw';
            else % Catalin's format (samples x channels)
                if ~isempty(obj.hCfg.loadTimeLimits)
                    nSamples = size(samplesRaw, 1);
                    lims = min(max(round(obj.hCfg.loadTimeLimits * obj.hCfg.sampleRate), 1), nSamples);

                    samplesRaw = samplesRaw(lims(1):lims(end), :);
                end

                if nargout > 1
                    channelMeans = single(mean(samplesRaw, 2)');
                end
            end
        end

        function [samplesOut, keepMe] = filterSamples(obj, samplesIn, windowPre, windowPost)
            %FILTERSAMPLES Denoise and filter raw samples, apply CAR
            tFilt = tic;

            fprintf('\tFiltering spikes...');
            samplesIn = [windowPre; samplesIn; windowPost];

            samplesIn_ = samplesIn; % keep a copy in CPU
            try
                samplesIn = jrclust.utils.tryGpuArray(samplesIn, obj.hCfg.useGPU);

                if obj.hCfg.fftThreshMAD > 0
                    samplesIn = jrclust.utils.fftClean(samplesIn, obj.hCfg.fftThreshMAD, obj.hCfg.ramToGPUFactor);
                end
            catch % GPU denoising failed, retry in CPU
                obj.hCfg.useGPU = false;

                samplesIn = samplesIn_;
                if obj.hCfg.fftThreshMAD > 0
                    samplesIn = jrclust.utils.fftClean(samplesIn, obj.hCfg.fftThreshMAD, obj.hCfg.ramToGPUFactor);
                end
            end

            % filter spikes; samples go in padded and come out padded
            try
                [samplesOut, channelMeans] = jrclust.utils.filtCar(samplesIn, [], [], false, obj.hCfg);
            catch % GPU filtering failed, retry in CPU
                obj.hCfg.useGPU = false;

                samplesIn = samplesIn_;
                [samplesOut, channelMeans] = jrclust.utils.filtCar(samplesIn, [], [], false, obj.hCfg);
            end

            clear rawSamples_;

            % common mode rejection
            if obj.hCfg.blank_thresh > 0
                if isempty(channelMeans) % carMode='whiten'
                    channelMeans = jrclust.utils.getCAR(samplesOut, obj.hCfg.carMode, obj.hCfg.ignoreSites);
                end

                keepMe = jrclust.utils.carReject(channelMeans(:), obj.hCfg.blank_period_ms, obj.hCfg.blank_thresh, obj.hCfg.sampleRate);
                fprintf('!! rejecting %0.3f %% of time due to motion !!', (1 - mean(keepMe))*100);
            else
                keepMe = [];
            end

            fprintf('\tdone (%0.2f) s\n', toc(tFilt));
        end

        function findSpikes(obj, samplesIn, keepMe, spTimes, spSites, siteThresh_ nPadPre, nPadPost)
            %FINDSPIKES detect spikes or use the one passed from the input (importing)
            if isempty(siteThresh_)
                try
                    siteThresh_ = jrclust.utils.tryGather(int16(mr2rms_(samplesIn, 1e5) * obj.hCfg.qqFactor));
                catch
                    obj.hCfg.fGpu = false;
                    siteThresh_ = int16(mr2rms_(jrclust.utils.tryGather(samplesIn), 1e5) * obj.hCfg.qqFactor);
                end
            end

            if isempty(spTimes) || isempty(spSites)
                try
                    obj.hCfg.nPad_pre = nPadPre;
                catch
                    obj.hCfg.addprop('nPad_pre');
                    obj.hCfg.nPad_pre = nPadPre;
                end
                [spTimes, spAmps, spSites] = jrclust.utils.detectSpikes(samplesIn, siteThresh_, keepMe, obj.hCfg);
            else
                spTimes = spTimes + nPadPre;
                spAmps = samplesIn(sub2ind(size(samplesIn), spTimes, spSites)); % @TODO read spikes at the site and time
            end

            % reject spikes within the overlap region
            if nPadPre > 0 || nPadPost > 0
                bounds = [nPadPre + 1, size(samplesIn, 1) - nPadPost]; % inclusive
                inBounds = spTimes >= bounds(1) & spTimes <= bounds(2);

                spTimes = spTimes(inBounds);
                spAmps = spAmps(inBounds);
                spSites = spSites(inBounds);
            end
            
            obj.siteThresh{end+1} = siteThresh_;
            obj.spikeTimes{end+1} = spTimes;
            obj.spikeAmps{end+1} = jrclust.utils.tryGather(spAmps);
            obj.spikeSites{end+1} = spSites;
        end

        function extractFeatures(obj, samplesRaw, samplesFilt, nPad_pre)
            % Extract spike waveforms and build a spike table
            fprintf('\tExtracting features');
            tf = tic;

            spTimes = obj.spikeTimes{end};
            spSites = obj.spikeSites{end};

            spSites_ = jrclust.utils.tryGpuArray(spSites);

            % spRaw, spFilt are nSamples x nChannels x nSpikes
            [spRaw, spFilt, spTimes] = obj.mn2tn_wav_(samplesRaw, samplesFilt, spSites_, spTimes);
            fprintf('.');

            if obj.hCfg.nFet_use >= 2
                spSites2_ = find_site_spk23_(spFilt, spSites_, obj.hCfg);
                spFilt2 = mn2tn_wav_spk2_(samplesFilt, spSites2_, spTimes, obj.hCfg);
            else
                [spSites2_, spFilt2] = deal([]);
            end

            % Cancel overlap
            if get_set_(obj.hCfg, 'fCancel_overlap', 0)
                try
                    [spFilt, spFilt2] = cancel_overlap_spk_(spFilt, spFilt2, spTimes, spSites, spSites2_, vnThresh_site, obj.hCfg);
                catch
                    fprintf(2, 'fCancel_overlap failed\n');
                end
            end

            obj.spikesRaw{end+1} = jrclust.utils.tryGather(spRaw);
            assert_(obj.hCfg.maxSite*2+1 - obj.hCfg.nSites_ref > 0, 'maxSite*2+1 - nSites_ref must be greater than 0');
            
            if obj.hCfg.nFet_use == 1
                mrFet1 = trWav2fet_(spFilt, obj.hCfg);
                fprintf('.');

                spFeatures = permute(mrFet1, [1,3,2]); %nSite x nFet x nSpk
                miSites = spSites_(:);
            elseif obj.hCfg.nFet_use == 2
                mrFet1 = trWav2fet_(spFilt, obj.hCfg);
                fprintf('.');

                mrFet2 = trWav2fet_(spFilt2, obj.hCfg);
                fprintf('.');

                spFeatures = permute(cat(3, mrFet1, mrFet2), [1, 3, 2]); % nSite x nFet x nSpk
                miSites = [spSites_(:), spSites2_(:)]; % nSpk x nFet
            else % obj.hCfg.nFet_use == 3
                [spSites2_, spSites3_] = find_site_spk23_(spFilt, spSites_, obj.hCfg);
                fprintf('.');

                mrFet1 = trWav2fet_(spFilt, obj.hCfg);
                fprintf('.');

                mrFet2 = trWav2fet_(spFilt2, obj.hCfg);
                fprintf('.');

                mrFet3 = trWav2fet_(mn2tn_wav_spk2_(samplesFilt, spSites3_, spTimes, obj.hCfg), obj.hCfg);
                fprintf('.');

                spFeatures = permute(cat(3, mrFet1, mrFet2, mrFet3), [1,3,2]); % nSite x nFet x nSpk
                miSites = [spSites_(:), spSites2_(:), spSites3_(:)]; % nSpk x nFet
            end

            if nPad_pre > 0
                spTimes = spTimes - nPad_pre;
            end

            obj.spikeTimes{end} = jrclust.utils.tryGather(spTimes);
            obj.spikesFilt{end+1} = jrclust.utils.tryGather(spFilt);
            obj.spikeFeatures{end+1} = jrclust.utils.tryGather(spFeatures);
            obj.miSite_spk{end+1} = jrclust.utils.tryGather(miSites);
            fprintf('done (%0.2f s)\n', toc(tf));
        end

        function [spikesRaw, spikesFilt, spTimes] = mn2tn_wav_(obj, samplesRaw, samplesFilt, spSites, spTimes)
            nSpks = numel(spTimes);
            nSites = numel(obj.hCfg.viSite2Chan);
            spkLim_wav = obj.hCfg.spkLim;
            spkLim_raw = obj.hCfg.spkLim_raw;
            nSites_spk = (obj.hCfg.maxSite * 2) + 1;

            % tensors, nSamples x nChannels x nSpikes
            spikesRaw = zeros(diff(spkLim_raw) + 1, nSites_spk, nSpks, 'like', samplesRaw);
            spikesFilt = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'like', samplesFilt);

            % Realignment parameters
            fRealign_spk = get_set_(obj.hCfg, 'fRealign_spk', 0); %0,1,2
            spTimes = jrclust.utils.tryGpuArray(spTimes, isGpu_(samplesRaw));
            spSites = jrclust.utils.tryGpuArray(spSites, isGpu_(samplesRaw));

            if isempty(spSites)
                spikesRaw = permute(mr2tr3_(samplesRaw, spkLim_raw, spTimes), [1,3,2]);
                spikesFilt = permute(mr2tr3_(samplesFilt, spkLim_wav, spTimes), [1,3,2]);
            else
                for iSite = 1:nSites
                    viiSpk11 = find(spSites == iSite);
                    if isempty(viiSpk11)
                        continue;
                    end

                    viTime_spk11 = spTimes(viiSpk11); % already sorted by time
                    viSite11 = obj.hCfg.siteNeighbors(:, iSite);

                    try
                        tnWav_spk1 = mr2tr3_(samplesFilt, spkLim_wav, viTime_spk11, viSite11);

                        if fRealign_spk == 1
                            [tnWav_spk1, viTime_spk11] = spkwav_realign_(tnWav_spk1, samplesFilt, spkLim_wav, viTime_spk11, viSite11, obj.hCfg);
                            spTimes(viiSpk11) = viTime_spk11;
                        elseif fRealign_spk == 2
                            tnWav_spk1 = spkwav_align_(tnWav_spk1, obj.hCfg);
                        end

                        spikesFilt(:, :, viiSpk11) = permute(tnWav_spk1, [1, 3, 2]);
                        spikesRaw(:,:,viiSpk11) = permute(mr2tr3_(samplesRaw, spkLim_raw, viTime_spk11, viSite11), [1,3,2]); %raw
                    catch % GPU failure
                        obj.isError = true;
                    end
                end
            end
        end
    end
end

