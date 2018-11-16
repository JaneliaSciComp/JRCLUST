function samplesOut = fftdiffFilter(samplesIn, freqLim, ramToGPUFactor)
    %FFTDIFFFILTER Summary of this function goes here
    if isempty(samplesIn)
        return;
    end

    useGPU = isa(samplesIn, 'gpuArray');

    nSamples = size(samplesIn, 1);
    [nLoads, nSamplesLoad, nSamplesLast] = jrclust.utils.partitionLoad(nSamples, round(nSamples/ramToGPUFactor));
    samplesOut = zeros(size(samplesIn), 'like', samplesIn);
    for iLoad = 1:nLoads
        offset = (iLoad-1) * nSamplesLoad;

        if iLoad < nLoads
            rows = offset + (1:nSamplesLoad);
        else
            rows = offset + (1:nSamplesLast);
        end
        samplesIn_ = samplesIn(rows, :);

        if useGPU
            try
                samplesOut(rows, :) = fftdOP(samplesIn_, freqLim);
            catch
                useGPU = 0;
            end
        end
        if ~useGPU % use CPU
            samplesOut(rows, :) = fftdOP(jrclust.utils.tryGather(samplesIn_), freqLim);
        end
    end % for
end

%% LOCAL FUNCTIONS
function samplesOut = fftdOP(samplesIn, freqLim)
    %FFTDOP Core of fftdiffFilter
    n = size(samplesIn, 1);

    n1 = round(n/2*freqLim(2));
    npow2 = 2^nextpow2(n);
    w = single(pi*1i)*single([linspace(0, 1, n1), linspace(1, -1, npow2-2*n1), linspace(-1, 0, n1)]');
    samplesOut = real(ifft(bsxfun(@times, fft(single(samplesIn), npow2), w), 'symmetric'));
    samplesOut = cast(samplesOut(1:n,:), class_(samplesIn));
end
