function [samplesOut, useGPU] = fftClean(samplesIn, fftThresh, ramToGPUFactor)
    %FFTCLEAN Remove high-frequency noise from samples
    if fftThresh == 0 || isempty(samplesIn)
        samplesOut = samplesIn;
        return;
    end

    useGPU = isa(samplesIn, 'gpuArray');

    fprintf('Applying FFT cleanup...');
    t = tic;

    nSamples = size(samplesIn, 1); % total number of samples (rows) in array
    [nLoads, nSamplesLoad, nSamplesFinal] = jrclust.utils.partitionLoad(nSamples, round(nSamples/ramToGPUFactor));
    samplesOut = zeros(size(samplesIn), 'like', samplesIn);

    for iLoad = 1:nLoads
        offset = (iLoad - 1)*nSamplesLoad;

        if iLoad < nLoads
            rows = offset + (1:nSamplesLoad);
        else
            rows = offset + (1:nSamplesFinal);
        end
        samplesOut_ = samplesIn(rows, :);

        if useGPU
            try
                samplesOut(rows, :) = doFFTClean(samplesOut_, fftThresh);
            catch
                useGPU = 0;
            end
        end

        if ~useGPU
            samplesOut(rows, :) = doFFTClean(jrclust.utils.tryGather(samplesOut_), fftThresh);
        end
    end % for

    fprintf(' done (%0.2f s)', toc(t));
end

%% LOCAL FUNCTIONS
function samplesOut = doFFTClean(samplesIn, fftThresh)
    %DOFFTCLEAN Actually perform the FFT clean
    nBins = 20;
    nSkipMed = 4; % skip this many samples when estimating median
    nw = 3;       % frequency neighbors to set to zero

    if fftThresh == 0 || isempty(fftThresh)
        samplesOut = samplesIn;
        return;
    end

    % get the next-largest power of 2 after nSamples (for padding)
    nSamples = size(samplesIn, 1);
    nSamplesPad = 2^nextpow2(nSamples);

    cls = class(samplesIn);

    for iRetry = 1:2
        try
            samplesOut = single(samplesIn);
            sampleMeans = mean(samplesOut, 1);
            samplesOut = bsxfun(@minus, samplesOut, sampleMeans); % center the samples

            if nSamples < nSamplesPad
                samplesOut = fft(samplesOut, nSamplesPad);
            else
                samplesOut = fft(samplesOut);
            end

            break; % success
        catch % failure, try again
            fprintf('!! GPU processing failed, retrying on CPU !!');
             samplesIn = jrclust.utils.tryGather(samplesIn);
        end
    end % for

    % find frequency outliers
    n1 = nSamplesPad/2; % previous power of 2
    viFreq = (1:n1)';
    vrFft = mean(bsxfun(@times, abs(samplesOut(1+viFreq,:)), viFreq), 2);
    n2 = round(n1/nBins);

    for iBin = 1:nBins
        vi1 = (n2*(iBin - 1):n2*iBin) + 1;
        if iBin == nBins
            vi1(vi1 > n1) = [];
        end

        % MAD transform
        vrFftMad = vrFft(vi1);
        vrFftMad = vrFftMad - median(vrFftMad(1:nSkipMed:end));
        vrFft(vi1) = vrFftMad / median(abs(vrFftMad(1:nSkipMed:end)));
    end

    % broaden spectrum
    isNoise = vrFft > fftThresh;
    vi_noise = find(isNoise);
    for i_nw = 1:nw
        viA = vi_noise - i_nw;
        viA(viA < 1)=[];

        viB = vi_noise + i_nw;
        viB(viB > n1)=[];

        isNoise(viA) = 1;
        isNoise(viB) = 1;
    end

    vi_noise = find(isNoise);
    samplesOut(1 + vi_noise, :) = 0;
    samplesOut(end - vi_noise + 1, :) = 0;

    % inverse transform back to the time domain
    samplesOut = real(ifft(samplesOut, nSamplesPad, 'symmetric'));

    if nSamples < nSamplesPad
        samplesOut = samplesOut(1:nSamples,:);
    end

    samplesOut = bsxfun(@plus, samplesOut, sampleMeans); % add mean back in
    samplesOut = cast(samplesOut, cls); % cast back to the original type
end
