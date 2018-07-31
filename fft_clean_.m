%--------------------------------------------------------------------------
% 2017/12/1 JJJ: auto-cast and memory division
function [mnWav1, useGPU] = fft_clean_(mnWav, P)
    useGPU = get_set_(P, 'useGPU', isGpu_(mnWav));
    if isempty(P.fft_thresh) || P.fft_thresh==0 || isempty(mnWav), mnWav1=mnWav; return; end
    [vcClass, useGPU_mnWav] = class_(mnWav);
    fprintf('Applying FFT cleanup...'); t1=tic;
    if 0
        if useGPU
            mnWav1 = fft_clean__(mnWav, P.fft_thresh);
        else
            mnWav1 = fft_clean__(gather_(mnWav), P.fft_thresh);
        end
    else
        nLoads_gpu = get_set_(P, 'nLoads_gpu', 8); % GPU load limit
        nSamples = size(mnWav,1);
        [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples, round(nSamples/nLoads_gpu));
        mnWav1 = zeros(size(mnWav), 'like', mnWav);
        for iLoad = 1:nLoad1
            iOffset = (iLoad-1) * nSamples_load1;
            if iLoad<nLoad1
                vi1 = (1:nSamples_load1) + iOffset;
            else
                vi1 = (1:nSamples_last1) + iOffset;
            end
            mnWav1_ = mnWav(vi1,:);
            if useGPU % use GPU
                try
                    if ~useGPU_mnWav, mnWav1_ = gpuArray_(mnWav1_); end
                    mnWav1(vi1,:) = fft_clean__(mnWav1_, P.fft_thresh);
                catch
                    useGPU = 0;
                end
            end
            if ~useGPU % use CPU
                mnWav1(vi1,:) = fft_clean__(gather_(mnWav1_), P.fft_thresh);
            end
        end %for
    end
    fprintf(' took %0.1fs', toc(t1));
end %func
