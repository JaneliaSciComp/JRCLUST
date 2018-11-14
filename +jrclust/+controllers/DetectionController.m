classdef DetectionController
    %DETECTIONCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        hCfg;
        hRecs;
        nRecs;
    end
    
    % LIFECYCLE
    methods
        function obj = DetectionController(hCfg)
            obj.hCfg = hCfg;
            obj.nRecs = numel(obj.hCfg.rawRecordings);
            obj.hRecs = jrclust.models.Recording.empty;
        end
    end

    % USER METHODS
    methods
        function res = detect(obj, givenTimes, givenSites)
            res = struct();
            t0 = tic;

            if nargin < 2
                givenTimes = [];
            end
            if nargin < 3
                givenSites = [];
            end

            givenTimes = givenTimes(:);
            givenSites = givenSites(:);

            dtype = obj.hCfg.dtype;
            nChans = obj.hCfg.nChans;
            headerOffset = obj.hCfg.headerOffset;

            % get manually set spike thresholds
            if ~isempty(obj.hCfg.threshFile)
                try
                    S = load(obj.hCfg.threshFile);
                    obj.siteThresh = S.vnThresh_site;
                    fprintf('Loaded %s\n', obj.hCfg.threshFile);
                catch
                    warning('Could not load threshFile %s', obj.hCfg.threshFile);
                    obj.siteThresh = [];
                end
            else
                obj.siteThresh = [];
            end

            % load from files
            for iRec = 1:obj.nRecs
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

                    fprintf('Processing %d/%d of file %d/%d...\n', iLoad, nLoads, iRec, obj.nRecs);
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
                    if ~isempty(givenTimes)
                        inInterval = (givenTimes >= sampOffset & givenTimes < sampOffset + nSamples);
                        intTimes = givenTimes(inInterval) - sampOffset + 1; % shift spike timing

                        % take sites assoc with times between limits
                        if ~isempty(givenSites)
                            intSites = givenSites(inInterval);
                        end
                    end

                    % increment sample offset
                    sampOffset = sampOffset + nSamples;

                    % load next N samples to ensure we don't miss any spikes at the boundary
                    if iLoad < nLoads && obj.hCfg.nSamplesPad > 0
                        windowPost = rec.readROI(obj.hCfg.siteMap, sampOffset:sampOffset+obj.hCfg.nSamplesPad-1);
                        windowPost = obj.regularize(windowPost);
                    else
                        windowPost = [];
                    end

                    %%%%%%%%%%%%%%

                    [tnWav_raw_, tnWav_spk_, trFet_spk_, miSite_spk{end+1}, viTime_spk{end+1}, vrAmp_spk{end+1}, siteThresh{end+1}, obj.hCfg.fGpu] ...
                        = wav2spk_(samplesRaw, channelMeans, intTimes, intSites, windowPre, windowPost);
                    write_spk_(tnWav_raw_, tnWav_spk_, trFet_spk_);
                    viTime_spk{end} = viTime_spk{end} + nSamples1;
                    nSamples1 = nSamples1 + nSamples;

                    if iLoad < nLoads
                        windowPre = samplesRaw(end-obj.hCfg.nSamplesPad+1:end, :);
                    end
                    clear mnWav11 vrWav_mean11;
                    nLoads = nLoads + 1;
                end
            end

            res.runtime = toc(t0);
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

        function samplesOut = filterSamples(obj, samplesIn, windowPre, windowPost)
            %FILTERSAMPLES Denoise and filter raw samples
            tFilt = tic;

            fprintf('\tFiltering spikes...');
            samplesIn = [windowPre; samplesIn; windowPost];

            samplesIn_ = samplesIn; % keep a copy in CPU
            try
                samplesIn = jrclust.utils.tryGpuArray(samplesIn, obj.hCfg.useGPU);

                if obj.hCfg.fft_thresh > 0
                    samplesIn = jrclust.utils.fftClean(samplesIn, P);
                end

                [samplesOut, vnWav11] = jrclust.utils.filtCar(samplesIn, P);
            catch % GPU failure
                obj.hCfg.useGPU = false;
                samplesIn = samplesIn_;
                if obj.hCfg.fft_thresh > 0
                    samplesIn = jrclust.utils.fftClean(samplesIn, P);
                end

                [samplesOut, vnWav11] = jrclust.utils.filtCar(samplesIn, P);
            end

            clear rawSamples_;

            % common mode rejection
            if obj.hCfg.blank_thresh > 0
                if isempty(vnWav11)
                    vnWav11 = mr2ref_(samplesOut, obj.hCfg.vcCommonRef, obj.hCfg.viSiteZero); %channelMeans(:);
                end

                vlKeep_ref = car_reject_(vnWav11(:), P);
                fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
            else
                vlKeep_ref = [];
            end

            % UNDO PADDING HERE?

            fprintf('\ttook %0.1fs\n', toc(tFilt));
        end

        %function [tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk, viTime_spk, vnAmp_spk, siteThresh, fGpu] = ...
        function w2sS = wav2spk_(obj, rawSamples, channelMeans, spikeTimes, spikeSites, windowPre, windowPost)
            %WAV2SPK_ wave to spike I guess
            [tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk] = deal([]);
            
            nSite_use = obj.hCfg.nSiteDir*2 - obj.hCfg.nSitesExcl + 1;
            if nSite_use == 1
                nFet_use = 1;
            else
                nFet_use = obj.hCfg.nFet_use;
            end

            % these don't appear to be used anywhere, may add back in later
            % fMerge_spk = 1;
            % fShift_pos = 0; % shift center position based on center of mass

            
            nSamplesPre = size(windowPre, 1);
            

            switch get_set_(P, 'vcFilter_detect', '')
                case {'', 'none'}, mnWav3 = mnWav2;
                case 'ndist'
                [mnWav3, nShift_post] = filter_detect_(rawSamples, P); % pass raw trace
                otherwise
                [mnWav3, nShift_post] = filter_detect_(mnWav2, P); % pass filtered trace
            end

            %-----
            % detect spikes or use the one passed from the input (importing)
            if isempty(obj.siteThresh)
                try
                    obj.siteThresh = jrclust.utils.tryGather(int16(mr2rms_(mnWav3, 1e5) * obj.hCfg.qqFactor));
                catch
                    obj.siteThresh = int16(mr2rms_(jrclust.utils.tryGather(mnWav3), 1e5) * obj.hCfg.qqFactor);
                    obj.hCfg.fGpu = 0;
                end
            end
            if isempty(spikeTimes) || isempty(spikeSites)
                P_ = setfield(P, 'nPad_pre', nSamplesPre);
                [spikeTimes, vnAmp_spk, spikeSites] = detect_spikes_(mnWav3, obj.siteThresh, vlKeep_ref, P_);
            else
                spikeTimes = spikeTimes + nSamplesPre;
                vnAmp_spk = mnWav3(sub2ind(size(mnWav3), spikeTimes, spikeSites)); % @TODO read spikes at the site and time
            end
            vnAmp_spk = jrclust.utils.tryGather(vnAmp_spk);
            % if nShift_post~=0, viTime_spk = viTime_spk + nShift_post; end % apply possible shift due to filtering

            % reject spikes within the overlap region
            if ~isempty(windowPre) || ~isempty(windowPost)
                ilim_spk = [nSamplesPre+1, size(mnWav3,1) - size(windowPost,1)]; %inclusive
                viKeep_spk = find(spikeTimes >= ilim_spk(1) & spikeTimes <= ilim_spk(2));
                [spikeTimes, vnAmp_spk, spikeSites] = multifun_(@(x)x(viKeep_spk), spikeTimes, vnAmp_spk, spikeSites);
            end %if
            if isempty(spikeTimes), return; end


            %-----
            % Extract spike waveforms and build a spike table
            fprintf('\tExtracting features'); t_fet = tic;
            % mnWav2 = jrclust.utils.tryGather(mnWav2); %do in CPU. 10.2s in GPU, 10.4s in CPU
            % if fRecenter_spk % center site is where the energy is the highest, if disabled min is chosen
            %     tnWav_spk = mn2tn_wav_spk2_(mnWav2, viSite_spk, viTime_spk, P);
            %     %[~, viMaxSite_spk] = max(squeeze_(std(single(tnWav_spk))));
            %     [~, viMaxSite_spk] = max(squeeze_(max(tnWav_spk) - min(tnWav_spk)));
            %     viSite_spk = obj.hCfg.miSites(sub2ind(size(obj.hCfg.miSites), viMaxSite_spk(:), viSite_spk));
            % end
            viSite_spk_ = jrclust.utils.tryGpuArray(spikeSites);
            [tnWav_spk_raw, tnWav_spk, spikeTimes] = mn2tn_wav_(rawSamples, mnWav2, viSite_spk_, spikeTimes, P); fprintf('.');
            if nFet_use >= 2
                viSite2_spk = find_site_spk23_(tnWav_spk, viSite_spk_, P);
                tnWav_spk2 = mn2tn_wav_spk2_(mnWav2, viSite2_spk, spikeTimes, P);
            else
                [viSite2_spk, tnWav_spk2] = deal([]);
            end

            %-----
            % Cancel overlap
            if get_set_(P, 'fCancel_overlap', 0)
                try
                    [tnWav_spk, tnWav_spk2] = cancel_overlap_spk_(tnWav_spk, tnWav_spk2, spikeTimes, spikeSites, viSite2_spk, obj.siteThresh, P);
                catch
                    fprintf(2, 'fCancel_overlap failed\n');
                end
            end

            tnWav_spk_raw = jrclust.utils.tryGather(tnWav_spk_raw);
            assert_(nSite_use >0, 'nSites_use = maxSite*2+1 - nSites_ref must be greater than 0');
            switch nFet_use
                case 3
                [viSite2_spk, viSite3_spk] = find_site_spk23_(tnWav_spk, viSite_spk_, P); fprintf('.');
                mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
                mrFet2 = trWav2fet_(tnWav_spk2, P); fprintf('.');
                mrFet3 = trWav2fet_(mn2tn_wav_spk2_(mnWav2, viSite3_spk, spikeTimes, P), P); fprintf('.');
                trFet_spk = permute(cat(3, mrFet1, mrFet2, mrFet3), [1,3,2]); %nSite x nFet x nSpk
                miSite_spk = [viSite_spk_(:), viSite2_spk(:), viSite3_spk(:)]; %nSpk x nFet
                case 2
                mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
                mrFet2 = trWav2fet_(tnWav_spk2, P); fprintf('.');
                trFet_spk = permute(cat(3, mrFet1, mrFet2), [1,3,2]); %nSite x nFet x nSpk
                miSite_spk = [viSite_spk_(:), viSite2_spk(:)]; %nSpk x nFet
                case 1
                mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
                trFet_spk = permute(mrFet1, [1,3,2]); %nSite x nFet x nSpk
                miSite_spk = [viSite_spk_(:)];
                otherwise
                error('wav2spk_: nFet_use must be 1, 2 or 3');
            end

            if nSamplesPre > 0, spikeTimes = spikeTimes - nSamplesPre; end
            [spikeTimes, trFet_spk, miSite_spk, tnWav_spk] = ...
            jrclust.utils.tryGather(spikeTimes, trFet_spk, miSite_spk, tnWav_spk);
            fGpu = obj.hCfg.fGpu;
            fprintf('\ttook %0.1fs\n', toc(t_fet));
        end %func
    end
end

