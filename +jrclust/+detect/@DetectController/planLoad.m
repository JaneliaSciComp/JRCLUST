function [nLoads, nSamplesLoad, nSamplesFinal] = planLoad(obj, hRec)
    %PLANLOAD Get number of samples to load in each chunk of a file
    nBytesSubset = subsetBytes(hRec, obj.hCfg.loadTimeLimits);

    bps = jrclust.utils.typeBytes(obj.hCfg.dataType);

    % nColumns in data matrix
    nSamples = floor(nBytesSubset / bps / obj.hCfg.nChans);

    % if not constrained by user, try to compute maximum bytes/load
    if isempty(obj.hCfg.maxBytesLoad)
        if obj.hCfg.useGPU
            S = gpuDevice(); % select first GPU device
            nBytes = floor(S(1).TotalMemory/2); % take half of total memory
        elseif ispc()
            S = memory();
            nBytes = floor(S.MaxPossibleArrayBytes);
        else % no hints given, assume 8 GiB
            nBytes = 2^33;
        end

        obj.hCfg.maxBytesLoad = floor(nBytes / obj.hCfg.gpuLoadFactor);
    end

    % if not constrained by user, try to compute maximum samples/load
    if isempty(obj.hCfg.maxSecLoad)
        nSamplesMax = floor(obj.hCfg.maxBytesLoad / obj.hCfg.nChans / bps);
    else
        nSamplesMax = floor(obj.hCfg.sampleRate * obj.hCfg.maxSecLoad);
    end

    if ~obj.hCfg.tallSkinny % load entire file, Catalin's format
        [nLoads, nSamplesLoad, nSamplesFinal] = deal(1, nSamples, nSamples);
    else
        [nLoads, nSamplesLoad, nSamplesFinal] = jrclust.utils.partitionLoad(nSamples, nSamplesMax);
    end
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