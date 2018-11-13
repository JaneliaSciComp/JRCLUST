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

            for iRec = 1:obj.nRecs
                fn = obj.hCfg.rawRecordings{iRec};
                rec = jrclust.models.Recording(fn, dtype, nChans, headerOffset);

                % get the number of bytes to load from recording
                nBytesFile = rec.subsetBytes(obj.hCfg.loadTimeLimits);
                % divide recording into many loads, `samples` are columns in data matrix
                [nLoads, nSamplesLoad, nSamplesFinal] = obj.planLoad(nBytesFile);

                % current column in data matrix
                sampOffset = 1;

                mnWav11_pre = [];
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
                    sampOffset = sampOffset + nSamples;

                    [mnWav11, vrWav_mean11] = obj.load_file_(samplesRaw);

                    fprintf('took %0.1fs\n', toc(t1));
                    if iLoad < nLoads
                        mnWav11_post = load_file_preview_(fid1, obj.hCfg);
                    else
                        mnWav11_post = [];
                    end

                    [viTime_spk11, viSite_spk11] = filter_spikes_(viTime_spk0, viSite_spk0, nSamples1 + [1, nSamples]);

                    [tnWav_raw_, tnWav_spk_, trFet_spk_, miSite_spk{end+1}, viTime_spk{end+1}, vrAmp_spk{end+1}, vnThresh_site{end+1}, obj.hCfg.fGpu] ...
                        = wav2spk_(mnWav11, vrWav_mean11, obj.hCfg, viTime_spk11, viSite_spk11, mnWav11_pre, mnWav11_post);
                    write_spk_(tnWav_raw_, tnWav_spk_, trFet_spk_);
                    viTime_spk{end} = viTime_spk{end} + nSamples1;
                    nSamples1 = nSamples1 + nSamples;
                    if iLoad < nLoads, mnWav11_pre = mnWav11(end-obj.hCfg.nPad_filt+1:end, :); end
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

        function [mnWav1, vrWav_mean1, dimm_wav] = load_file_(obj, samplesRaw)
            fSingle = 0; % output single
            if obj.hCfg.fTranspose_bin
                dimm_wav = [obj.hCfg.nChans, nSamples_load1];
            else % Catalin's format
                dimm_wav = [nSamples_load1, obj.hCfg.nChans];
            end
%             mnWav1 = fread_(fid_bin, dimm_wav, obj.hCfg.dtype);
            switch(obj.hCfg.dtype)
                case 'uint16'
                    mnWav1 = int16(single(samplesRaw)-2^15);

                case {'single', 'double'}
                    mnWav1 = int16(samplesRaw / obj.hCfg.bitScaling);
            end

            % flip the polarity
            if obj.hCfg.fInverse_file
                mnWav1 = -mnWav1;
            end

            % extract channels
            if obj.hCfg.fTranspose_bin
                vrWav_mean1 = single(mean(mnWav1, 1)); %6x faster to transpose in dimm1
                mnWav1 = mnWav1';
            else %Catalin's format. time x nChans
                if ~isempty(obj.hCfg.viSite2Chan), mnWav1 = mnWav1(:,obj.hCfg.viSite2Chan); end
                if ~isempty(obj.hCfg.tlim_load)
                    nSamples = size(mnWav1,1);
                    nlim_load = min(max(round(obj.hCfg.tlim_load * obj.hCfg.sRateHz), 1), nSamples);
                    mnWav1 = mnWav1(nlim_load(1):nlim_load(end), :);
                end
                vrWav_mean1 = single(mean(mnWav1, 2)');
                if fSingle
                    mnWav1 = single(mnWav1) * obj.hCfg.bitScaling;
                end
            end
            if ~isempty(vcFile)
                fclose(fid_bin);
            end
        end %func


    end
end

