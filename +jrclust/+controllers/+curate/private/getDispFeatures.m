function dispFeatures = getDispFeatures(hClust, hCfg, sitesToShow, selected, dispFeature)
    %GETDISPFEATURES Get features to show in feature projection view
    %   input:  sitesToShow, vector of site indices
    %   input:  selected, vector of primary (and optionally secondary)
    %           selected clusters
    %   input:  dispFeature, optional string, the feature to display
    %   output: dispFeatures, struct containing x-y values of background,
    %           foreground, and secondary foreground features
    if nargin < 5
        dispFeature = hCfg.dispFeature;
    end

    iCluster = selected(1);
    if numel(selected) == 2
        jCluster = selected(2);
    else
        jCluster = [];
    end

    % select subset of spikes
    spikesToShow = find(ismember(hClust.spikeSites, sitesToShow));
    timesToShow = hClust.spikeTimes(spikesToShow);

    % if we're only interested in a subset of time, just plot those spikes
    if ~isempty(hCfg.tLimFigProj)
        timeBound = round(hCfg.tLimFigProj*hCfg.sampleRate);
        inLimit = (timesToShow >= timeBound(1) & timesToShow <= timeBound(end));
        spikesToShow = spikesToShow(inLimit);
    end

    % get cluster assignments of spikes to show
    spikeClustersShow = hClust.spikeClusters(spikesToShow);

    bgSpikes = jrclust.utils.subsample(spikesToShow, 2*hCfg.nShow_proj);
    fgSpikes = jrclust.utils.subsample(spikesToShow(spikeClustersShow == iCluster), hCfg.nShow_proj);
    if ~isempty(jCluster)
        fgSpikes2 = randomSelect_(spikesToShow(spikeClustersShow == jCluster), hCfg.nShow_proj);
    else
        [fg2Y, fg2X] = deal([]);
    end

    if strcmp(dispFeature, 'pca')
        % compute first 2 principal vectors (on each of sitesToShow) from
        % all spikes in iCluster
        sampledWindows = permute(jrclust.utils.getSampledWindows(hClust, hClust.spikesByCluster{iCluster}, sitesToShow, false), [1, 3, 2]);
        [prVecs1, prVecs2] = jrclust.features.getPVSpikes(sampledWindows);

        % project background and foreground spikes onto principal vectors
        bgWindows = permute(jrclust.utils.getSampledWindows(hClust, bgSpikes, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
        [bgY, bgX] = jrclust.features.pcProjectSpikes(bgWindows, prVecs1, prVecs2);

        fgWindows = permute(jrclust.utils.getSampledWindows(hClust, fgSpikes, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
        [fgY, fgX] = jrclust.features.pcProjectSpikes(fgWindows, prVecs1, prVecs2);

        if ~isempty(jCluster)
            fgWindows2 = permute(jrclust.utils.getSampledWindows(hClust, fgSpikes2, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
            [fg2Y, fg2X] = jrclust.features.pcProjectSpikes(fgWindows2, prVecs1, prVecs2);
        end
    elseif strcmp(dispFeature, 'ppca')
        % compute first principal vectors (on each of sitesToShow) from
        % spikes in iCluster and jCluster, respectively
        [prVecs1, prVecs2] = jrclust.features.getPVClusters(sitesToShow, iCluster, jCluster);

        % project background and foreground spikes onto principal vectors
        bgWindows = permute(jrclust.utils.getSampledWindows(hClust, bgSpikes, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
        [bgY, bgX] = jrclust.features.pcProjectSpikes(bgWindows, prVecs1, prVecs2);

        fgWindows = permute(jrclust.utils.getSampledWindows(hClust, fgSpikes, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
        [fgY, fgX] = jrclust.features.pcProjectSpikes(fgWindows, prVecs1, prVecs2);

        if ~isempty(jCluster)
            fgWindows2 = permute(jrclust.utils.getSampledWindows(hClust, fgSpikes2, sitesToShow, false), [1, 3, 2]); % nSamples x nSpikes x nSites
            [fg2Y, fg2X] = jrclust.features.pcProjectSpikes(fgWindows2, prVecs1, prVecs2);
        end
    elseif strcmp(dispFeature, 'cov')
        [bgY, bgX] = getSpikeCov(hClust, bgSpikes, sitesToShow);
        [fgY, fgX] = getSpikeCov(hClust, fgSpikes, sitesToShow);

        if ~isempty(jCluster)
            [fg2Y, fg2X] = getSpikeCov(hClust, fgSpikes2, sitesToShow);
        end
    elseif strcmp(dispFeature, 'vpp')
        bgWindows = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust, bgSpikes, sitesToShow, false), hCfg);
        bgY = abs(permute(min(bgWindows), [2, 3, 1]));
        bgX = abs(permute(max(bgWindows), [2, 3, 1]));

        fgWindows = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust, fgSpikes, sitesToShow, false), hCfg);
        fgY = abs(permute(min(fgWindows), [2, 3, 1]));
        fgX = abs(permute(max(fgWindows), [2, 3, 1]));

        if ~isempty(jCluster)
            fgWindows2 = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust, fgSpikes2, sitesToShow, false), hCfg);
            fg2Y = abs(permute(min(fgWindows2), [2, 3, 1]));
            fg2X = abs(permute(max(fgWindows2), [2, 3, 1]));
        end
    end

    dispFeatures = struct('bgY', abs(bgY), ...
                          'bgX', abs(bgX), ...
                          'fgY', abs(fgY), ...
                          'fgX', abs(fgX), ...
                          'fg2Y', abs(fg2Y), ...
                          'fg2X', abs(fg2X));
end