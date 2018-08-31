%--------------------------------------------------------------------------
% 17/12/1 JJJ: Load size is not limited by FFT cleanup process (fft_thresh>0)
function [nLoad1, nSamples_load1, nSamples_last1] = planLoad(nBytes_file, P)
    % plan load file size according to the available memory and file size (nBytes_file)
    LOAD_FACTOR = 5; %GPU memory usage factor. 4x means 1/4 of GPU memory can be loaded

    nSamples1 = floor(nBytes_file/(bytesPerSample_(P.dataType)*P.nChans));

    if ~isfield(P, 'MAX_BYTES_LOAD') || isempty(P.MAX_BYTES_LOAD)
        try
            if P.useGPU
                S = gpuDevice(); % does not reset GPU
                memMax = floor(S(1).AvailableMemory());
            else
                S = memory();
                memMax = floor(S.MaxPossibleArrayBytes());
            end
        catch
            memMax = inf; % assume infinite memory
        end
        P.MAX_BYTES_LOAD = floor(memMax/LOAD_FACTOR);
    end

    if ~isfield(P, 'MAX_LOAD_SEC') || isempty(P.MAX_LOAD_SEC)
        nSamples_max = floor(P.MAX_BYTES_LOAD/(P.nChans*bytesPerSample_(P.dataType)));
    else
        nSamples_max = floor(P.sampleRateHz * P.MAX_LOAD_SEC);
    end

    if ~P.fTranspose_bin % load all in one, Catalin's format
        [nLoad1, nSamples_load1, nSamples_last1] = deal(1, nSamples1, nSamples1);
    else
        [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max);
    end
end % function
