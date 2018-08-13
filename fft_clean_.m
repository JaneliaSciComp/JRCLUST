%--------------------------------------------------------------------------
% 2017/12/1 JJJ: auto-cast and memory division
function [cleanedSamplesMatrix, useGPU] = fft_clean_(rawTraces, P)
    useGPU = get_set_(P, 'useGPU', isGpu_(rawTraces));
    if ~isfield(P, 'fft_thresh') || isempty(P.fft_thresh) || P.fft_thresh==0 || isempty(rawTraces)
        cleanedSamplesMatrix = rawTraces;
        return;
    end

    [vcClass, useGPU_mnWav] = class_(rawTraces);
    fprintf('Applying FFT cleanup...');

    t1 = tic;
    nLoads_gpu = get_set_(P, 'nLoads_gpu', 8); % GPU load limit
    nSamples = size(rawTraces,1);
    [nLoads, nSamples_load1, nSamples_last1] = partition_load_(nSamples, round(nSamples/nLoads_gpu));

    cleanedSamplesMatrix = zeros(size(rawTraces), 'like', rawTraces);

    for iLoad = 1:nLoads
        iOffset = (iLoad - 1)*nSamples_load1;
        if iLoad<nLoads
            vi1 = (1:nSamples_load1) + iOffset;
        else
            vi1 = (1:nSamples_last1) + iOffset;
        end

        mnWav1_ = rawTraces(vi1,:);

        if useGPU % use GPU
            try
                if ~useGPU_mnWav, mnWav1_ = gpuArray_(mnWav1_); end
                cleanedSamplesMatrix(vi1,:) = fft_clean__(mnWav1_, P.fft_thresh);
            catch
                useGPU = 0;
            end
        end

        if ~useGPU % use CPU
            cleanedSamplesMatrix(vi1,:) = fft_clean__(gather_(mnWav1_), P.fft_thresh);
        end
    end %for
    fprintf(' took %0.1fs', toc(t1));
end %func
