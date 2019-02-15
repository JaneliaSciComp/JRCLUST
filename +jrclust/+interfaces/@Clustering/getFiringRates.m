function rates = getFiringRates(obj, clusters, nSamples)
    %GETFIRINGRATES Compute firing rates for specified clusters
    if nargin < 2
        clusters = 1:obj.nClusters;
    end
    if nargin < 3 || isempty(nSamples)
        nSamples = round(obj.hCfg.recDurationSec()*obj.hCfg.frSampleRate);
    end

    nFilt = round(obj.hCfg.frSampleRate * obj.hCfg.frPeriod / 2);
    if strcmp(obj.hCfg.frFilterShape, 'triangle')
        filtKernel = ([1:nFilt, nFilt-1:-1:1]'/nFilt*2/obj.hCfg.frPeriod);
    elseif strctmp(obj.hCfg.frFilterShape, 'rectangle')
        filtKernel = (ones(nFilt*2, 1) / obj.hCfg.frPeriod);
    end % switch

    filtKernel = single(filtKernel);

    rates = zeros([nSamples, numel(clusters)], 'single');
    for iiCluster = 1:numel(clusters)
        iCluster = clusters(iiCluster);
        clusterSpikes = obj.spikesByCluster{iCluster};
        clusterTimes = obj.spikeTimes(clusterSpikes);

        dsTimes = round(obj.hCfg.frSampleRate*double(clusterTimes)/obj.hCfg.sampleRate);
        dsTimes = max(min(dsTimes, nSamples), 1);

        rates(dsTimes, iiCluster) = 1;
        rates(:, iiCluster) = conv(rates(:, iiCluster), filtKernel, 'same');
    end
end