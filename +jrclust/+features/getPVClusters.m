function [prVecs1, prVecs2] = getPVClusters(hClust, sites, iCluster, jCluster)
    %GETPVCLUSTERS Compute principal vectors on sites from one or two clusters
    %   input:  sites, a vector of site indices
    %   input:  iCluster, primary cluster ID of interest
    %   input:  jCluster, optional, secondary cluster ID of interest
    %   output: prVecs1, nSamples x nSites; 1st principal vectors of
    %           iCluster spike waveforms on sites
    %   output: prVecs2, nSamples x nSites; if jCluster is given, 1st
    %           principal vectors of jCluster spike waveforms on sites;
    %           if jCluster is not given, 2nd principal vectors of iCluster
    %           spike waveforms on sites

    if nargin < 4
        jCluster = [];
    end

    MAX_SAMPLE = 10000;
    sampledSpikes = jrclust.utils.subsample(hClust.spikesByCluster{iCluster}, MAX_SAMPLE);
    sampledWindows = permute(hClust.getSpikeWindows(sampledSpikes, sites, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites

    nSamples = size(sampledWindows, 1);
    nSites = numel(sites);

    [prVecs1, prVecs2] = deal(zeros(nSamples, nSites, 'single'));
    if isempty(jCluster)
        % compute 1st and 2nd principal vectors off of the spikes from
        % iCluster and store them in prVecs{1,2} respectively
        for kSite = 1:nSites
            prVecs = jrclust.features.getPVSamples(sampledWindows(:, :, kSite));
            prVecs1(:, kSite) = prVecs(:, 1);
            prVecs2(:, kSite) = prVecs(:, 2);
        end
    else
        % compute 1st principal vectors off of the spikes from iCluster and
        % jCluster and store them in prVecs{1,2} respectively
        sampledSpikes2 = jrclust.utils.subsample(hClust.spikesByCluster{jCluster}, MAX_SAMPLE);
        sampledWindows2 = permute(hClust.getSpikeWindows(sampledSpikes2, sites, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
        for kSite = 1:nSites
            iVecs = jrclust.features.getPVSamples(sampledWindows(:, :, kSite));
            jVecs = jrclust.features.getPVSamples(sampledWindows2(:, :, kSite));
            prVecs1(:, kSite) = iVecs(:, 1);
            prVecs2(:, kSite) = jVecs(:, 1);
        end
    end
end
