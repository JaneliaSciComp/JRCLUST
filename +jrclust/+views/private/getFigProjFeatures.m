function dispFeatures = getFigProjFeatures(hClust, sitesToShow, selected)
    %GETFIGPROJFEATURES Get features to show in feature projection view
    %   input:  sitesToShow, vector of site indices
    %   input:  selected, vector of primary (and optionally secondary)
    %           selected clusters
    %   output: dispFeatures, struct containing x-y values of background,
    %           foreground, and secondary foreground features
    hCfg = hClust.hCfg;
    dispFeature = hCfg.dispFeature;
    if strcmp(dispFeature, 'template') && ~isa(hClust, 'jrclust.sort.TemplateClustering')
        dispFeature = 'vpp';
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
    if ~isempty(hCfg.projTimeLimits)
        timeBound = round(hCfg.projTimeLimits*hCfg.sampleRate);
        inLimit = (timesToShow >= timeBound(1) & timesToShow <= timeBound(end));
        spikesToShow = spikesToShow(inLimit);
    end

    % get cluster assignments of spikes to show
    spikeClustersShow = hClust.spikeClusters(spikesToShow);

    bgSpikes = jrclust.utils.subsample(spikesToShow, 2*hCfg.nSpikesFigProj);
    fgSpikes = jrclust.utils.subsample(spikesToShow(spikeClustersShow == iCluster), hCfg.nSpikesFigProj);
    if ~isempty(jCluster)
        fg2Spikes = jrclust.utils.subsample(spikesToShow(spikeClustersShow == jCluster), hCfg.nSpikesFigProj);
    else
        [fg2Spikes, fg2YData, fg2XData] = deal([]);
    end

    switch dispFeature
        case 'pca'
            % compute first 2 principal vectors (on each of sitesToShow) from
            % all spikes in iCluster
            sampledWindows = permute(hClust.getSpikeWindows(hClust.spikesByCluster{iCluster}, sitesToShow, 0, 0), [1, 3, 2]);
            if all(hCfg.pcPair == [1 2])
                [prVecs1, prVecs2] = jrclust.features.getPVSpikes(sampledWindows);
            elseif all(hCfg.pcPair == [1 3])
                [prVecs1, ~, prVecs2] = jrclust.features.getPVSpikes(sampledWindows);
            else % [2 3]
                [~, prVecs1, prVecs2] = jrclust.features.getPVSpikes(sampledWindows);
            end

            % project background and foreground spikes onto principal vectors
            bgWindows = permute(hClust.getSpikeWindows(bgSpikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [bgYData, bgXData] = jrclust.features.pcProjectSpikes(bgWindows, prVecs1, prVecs2);

            fgWindows = permute(hClust.getSpikeWindows(fgSpikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [fgYData, fgXData] = jrclust.features.pcProjectSpikes(fgWindows, prVecs1, prVecs2);

            if ~isempty(jCluster)
                fgWindows2 = permute(hClust.getSpikeWindows(fg2Spikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
                [fg2YData, fg2XData] = jrclust.features.pcProjectSpikes(fgWindows2, prVecs1, prVecs2);
            end

        case 'ppca'
            % compute first principal vectors (on each of sitesToShow) from
            % spikes in iCluster and jCluster, respectively
            [prVecs1, prVecs2] = jrclust.features.getPVClusters(hClust, sitesToShow, iCluster, jCluster);

            % project background and foreground spikes onto principal vectors
            bgWindows = permute(hClust.getSpikeWindows(bgSpikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [bgYData, bgXData] = jrclust.features.pcProjectSpikes(bgWindows, prVecs1, prVecs2);

            fgWindows = permute(hClust.getSpikeWindows(fgSpikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [fgYData, fgXData] = jrclust.features.pcProjectSpikes(fgWindows, prVecs1, prVecs2);

            if ~isempty(jCluster)
                fgWindows2 = permute(hClust.getSpikeWindows(fg2Spikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
                [fg2YData, fg2XData] = jrclust.features.pcProjectSpikes(fgWindows2, prVecs1, prVecs2);
            end

        case 'cov'
            [bgYData, bgXData] = getSpikeCov(hClust, bgSpikes, sitesToShow);
            [fgYData, fgXData] = getSpikeCov(hClust, fgSpikes, sitesToShow);

            if ~isempty(jCluster)
                [fg2YData, fg2XData] = getSpikeCov(hClust, fg2Spikes, sitesToShow);
            end

        case 'vpp'
            bgWindows = hClust.getSpikeWindows(bgSpikes, sitesToShow, 0, 1); % use voltages 
            bgYData = abs(permute(min(bgWindows), [2, 3, 1]));
            bgXData = abs(permute(max(bgWindows), [2, 3, 1]));

            fgWindows = hClust.getSpikeWindows(fgSpikes, sitesToShow, 0, 1); % use voltages 
            fgYData = abs(permute(min(fgWindows), [2, 3, 1]));
            fgXData = abs(permute(max(fgWindows), [2, 3, 1]));

            if ~isempty(jCluster)
                fgWindows2 = hClust.getSpikeWindows(fg2Spikes, sitesToShow, 0, 1); % use voltages
                fg2YData = abs(permute(min(fgWindows2), [2, 3, 1]));
                fg2XData = abs(permute(max(fgWindows2), [2, 3, 1]));
            end

        case 'template' % currently Kilosort only
            if all(hCfg.pcPair == [1 2])
                [bgYData, bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [fgYData, fgXData] = hClust.pcFeaturesBySpike(fgSpikes, sitesToShow);
                if ~isempty(jCluster)
                    [fg2YData, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            elseif all(hCfg.pcPair == [1 3])
                [bgYData, ~, bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [fgYData, ~, fgXData] = hClust.pcFeaturesBySpike(fgSpikes, sitesToShow);
                if ~isempty(jCluster)
                    [fg2YData, ~, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            else % [2 3]
                [~, bgYData,bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [~, fgYData, fgXData] = hClust.pcFeaturesBySpike(fgSpikes, sitesToShow);
                if ~isempty(jCluster)
                    [~, fg2YData, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            end
    end

    dispFeatures = struct('bgYData', bgYData, ...
                          'bgXData', bgXData, ...
                          'fgYData', fgYData, ...
                          'fgXData', fgXData, ...
                          'fgSpikes', fgSpikes, ... % for finding spikes to split off
                          'fg2YData', fg2YData, ...
                          'fg2XData', fg2XData, ...
                          'fg2Spikes', fg2Spikes);
end