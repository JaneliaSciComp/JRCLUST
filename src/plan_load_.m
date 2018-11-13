%--------------------------------------------------------------------------
% 17/12/1 JJJ: Load size is not limited by FFT cleanup process (fft_thresh>0)
function [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file, P)
    % plan load file size according to the available memory and file size (nBytes_file1)
    LOAD_FACTOR = 5; %GPU memory usage factor. 4x means 1/4 of GPU memory can be loaded

    nSamples1 = floor(nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans);
    % nSamples_max = floor(mem_max_(P) / P.nChans / 4); % Bound by MAX_BYTES_LOAD
    if ~isfield(P, 'MAX_BYTES_LOAD'), P.MAX_BYTES_LOAD = []; end
    if isempty(P.MAX_BYTES_LOAD), P.MAX_BYTES_LOAD = floor(mem_max_(P) / LOAD_FACTOR); end
    if isempty(P.MAX_LOAD_SEC)
        nSamples_max = floor(P.MAX_BYTES_LOAD / P.nChans / bytesPerSample_(P.vcDataType));
    else
        nSamples_max = floor(P.sRateHz * P.MAX_LOAD_SEC);
    end

    if ~P.fTranspose_bin %load all in one, Catalin's format
        [nLoad1, nSamples_load1, nSamples_last1] = deal(1, nSamples1, nSamples1);
    else
        [nLoad1, nSamples_load1, nSamples_last1] = jrclust.utils.partitionLoad(nSamples1, nSamples_max);
    end
end %func
