function sampleRMS = estimateRMS(samplesIn, nSubsamples)
    %ESTIMATERMS Estimate RMS of samples using median
    if nargin < 2 || isempty(nSubsamples)
        sampleRMS = median(abs(samplesIn));
    else
        sampleRMS = median(abs(jrclust.utils.subsample(samplesIn, nSubsamples, 1)));
    end

    sampleRMS = single(sampleRMS) / 0.6745;
end
