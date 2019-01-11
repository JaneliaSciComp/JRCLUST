function [nLoads, nSamplesLoad, nSamplesFinal] = planLoad(nBytesFile, hCfg)
    %PLANLOAD Get number of samples to load in each chunk of a file
    bps = jrclust.utils.typeBytes(hCfg.dtype);

    % nColumns in data matrix
    nSamples = floor(nBytesFile / bps / hCfg.nChans);

    % if not constrained by user, try to compute maximum bytes/load
    if isempty(hCfg.maxBytesLoad)
        if hCfg.useGPU
            S = gpuDevice(); % select first GPU device
            nBytes = floor(S(1).TotalMemory/2); % take 1/2 of total memory
        elseif ispc()
            S = memory();
            nBytes = floor(S.MaxPossibleArrayBytes);
        else % no hints given, assume infinite memory
            nBytes = inf;
        end

        hCfg.maxBytesLoad = floor(nBytes / hCfg.gpuLoadFactor);
    end

    % if not constrained by user, try to compute maximum samples/load
    if isempty(hCfg.maxSecLoad)
        nSamplesMax = floor(hCfg.maxBytesLoad / hCfg.nChans / bps);
    else
        nSamplesMax = floor(hCfg.sampleRate * hCfg.maxSecLoad);
    end

    if ~hCfg.fTranspose_bin % load entire file, Catalin's format
        [nLoads, nSamplesLoad, nSamplesFinal] = deal(1, nSamples, nSamples);
    else
        [nLoads, nSamplesLoad, nSamplesFinal] = jrclust.utils.partitionLoad(nSamples, nSamplesMax);
    end
end