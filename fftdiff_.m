%--------------------------------------------------------------------------
function mnWav1 = fftdiff_(mnWav, P)

    useGPU = isGpu_(mnWav);
    nLoads_gpu = get_set_(P, 'nLoads_gpu', 8); % GPU load limit

    % [useGPU, nLoads_gpu] = deal(0, 1); %debug

    nSamples = size(mnWav,1);
    [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples, round(nSamples/nLoads_gpu));
    mnWav1 = zeros(size(mnWav), 'like', mnWav);
    freqLim_ = P.freqLim / (P.sRateHz / 2);
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
                mnWav1(vi1,:) = fftdiff__(mnWav1_, freqLim_);
            catch
                useGPU = 0;
            end
        end
        if ~useGPU % use CPU
            mnWav1(vi1,:) = fftdiff__(gather_(mnWav1_), freqLim_);
        end
    end %for
end %func
