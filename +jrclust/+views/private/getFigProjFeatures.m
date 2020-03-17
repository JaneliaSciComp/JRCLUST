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

    siteMask = ismember(hClust.spikeSites, sitesToShow);
    % if we're only interested in a subset of time, just plot those spikes
    if ~isempty(hCfg.projTimeLimits)
        timeBound = round(hCfg.projTimeLimits*hCfg.sampleRate);
        timeMask = (timesToShow >= timeBound(1) & timesToShow <= timeBound(end));
    else
        timeMask = true(size(hClust.spikeClusters));
    end

    % spikes occurring on these sites and within these times and NOT in iCluster
    bgMask = (hClust.spikeClusters ~= iCluster) & siteMask & timeMask;

    if ~isempty(jCluster)
        % spikes occurring on these sites and within these times and in jCluster
        jMask = (hClust.spikeClusters == jCluster) & siteMask & timeMask;

        % update bgMask to exclude spikes from jCluster
        bgMask = bgMask & (~jMask);
    else
        jMask = false(size(hClust.spikeClusters));
    end

    % subselect spikes from jCluster and background
    bgSpikes = jrclust.utils.subsample(find(bgMask), 2*hCfg.nSpikesFigProj);
    if any(jMask)
        fg2Spikes = jrclust.utils.subsample(find(jMask), hCfg.nSpikesFigProj);
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

            % get features for ALL foreground spikes on sitesToShow
            fgWindows = permute(hClust.getSpikeWindows(hClust.spikesByCluster{iCluster}, ...
                                                       sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
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

            % get features for ALL foreground spikes on sitesToShow
            fgWindows = permute(hClust.getSpikeWindows(hClust.spikesByCluster{iCluster}, ...
                                                       sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [fgYData, fgXData] = jrclust.features.pcProjectSpikes(fgWindows, prVecs1, prVecs2);

            if ~isempty(jCluster)
                fgWindows2 = permute(hClust.getSpikeWindows(fg2Spikes, sitesToShow, 0, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
                [fg2YData, fg2XData] = jrclust.features.pcProjectSpikes(fgWindows2, prVecs1, prVecs2);
            end

        case 'cov'
            [bgYData, bgXData] = getSpikeCov(hClust, bgSpikes, sitesToShow);
            [fgYData, fgXData] = getSpikeCov(hClust, hClust.spikesByCluster{iCluster}, sitesToShow);

            if ~isempty(jCluster)
                [fg2YData, fg2XData] = getSpikeCov(hClust, fg2Spikes, sitesToShow);
            end

        case 'vpp'
            bgWindows = hClust.getSpikeWindows(bgSpikes, sitesToShow, 0, 1); % use voltages
            bgYData = abs(permute(min(bgWindows), [2, 3, 1]));
            bgXData = abs(permute(max(bgWindows), [2, 3, 1]));

            % get features for ALL foreground spikes on sitesToShow
            fgWindows = hClust.getSpikeWindows(hClust.spikesByCluster{iCluster}, ...
                                               sitesToShow, 0, 1); % use voltages
            fgYData = abs(permute(min(fgWindows), [2, 3, 1])); % nSitesToShow x nSpikes
            fgXData = abs(permute(max(fgWindows), [2, 3, 1]));

            if ~isempty(jCluster)
                fgWindows2 = hClust.getSpikeWindows(fg2Spikes, sitesToShow, 0, 1); % use voltages
                fg2YData = abs(permute(min(fgWindows2), [2, 3, 1]));
                fg2XData = abs(permute(max(fgWindows2), [2, 3, 1]));
            end

        case 'template' % currently Kilosort only
            if all(hCfg.pcPair == [1 2])
                [bgYData, bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [fgYData, fgXData] = hClust.pcFeaturesBySpike(hClust.spikesByCluster{iCluster}, sitesToShow);
                if ~isempty(jCluster)
                    [fg2YData, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            elseif all(hCfg.pcPair == [1 3])
                [bgYData, ~, bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [fgYData, ~, fgXData] = hClust.pcFeaturesBySpike(hClust.spikesByCluster{iCluster}, sitesToShow);
                if ~isempty(jCluster)
                    [fg2YData, ~, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            else % [2 3]
                [~, bgYData,bgXData] = hClust.pcFeaturesBySpike(bgSpikes, sitesToShow);
                [~, fgYData, fgXData] = hClust.pcFeaturesBySpike(hClust.spikesByCluster{iCluster}, sitesToShow);
                if ~isempty(jCluster)
                    [~, fg2YData, fg2XData] = hClust.pcFeaturesBySpike(fg2Spikes, sitesToShow);
                end
            end
    end

    dispFeatures = struct('bgYData', bgYData, ...
                          'bgXData', bgXData, ...
                          'fgYData', fgYData, ...
                          'fgXData', fgXData, ...
                          'fg2YData', fg2YData, ...
                          'fg2XData', fg2XData);
end
