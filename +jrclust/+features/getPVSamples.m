function prVecs = getPVSamples(samplesIn)
    %GETPVSAMPLES Get top 3 principal vectors for spike waveforms
    %   samplesIn is nSamples x nSpikes
    MAX_SAMPLE = 10000;

    % select at most MAX_SAMPLE spikes to estimate the covariance matrix
    samplesIn = jrclust.utils.subsample(samplesIn, MAX_SAMPLE, 2);
    samplesIn = jrclust.utils.meanSubtract(single(samplesIn));

    % compute eigenvectors of covariance matrix
    nSamples = size(samplesIn, 2);
    covMat = (samplesIn*samplesIn')/(nSamples-1);
    [eigVecs, ~] = eig(covMat);

    % eigenvectors corresponding to 3 largest eigenvalues
    prVecs = eigVecs(:, end:-1:end-2);
end
