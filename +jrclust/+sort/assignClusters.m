function sRes = assignClusters(dRes, sRes, hCfg)
    %ASSIGNCLUSTERS Given rho-delta information, assign spikes to clusters
    sRes = computeCenters(dRes, sRes, hCfg);
    sRes.spikeClusters = [];

    % do initial assignment
    hCfg.updateLog('assignClusters', sprintf('Assigning clusters (nClusters: %d)', numel(sRes.clusterCenters)), 1, 0);
    sRes = doAssign(dRes, sRes, hCfg);
    nClusters = numel(sRes.clusterCenters);

    % split clusters by site distance
    sRes = splitDist(dRes, sRes, hCfg);

    % reassign spikes if we had to split
    if nClusters ~= numel(sRes.clusterCenters)
        sRes.spikeClusters = [];
        sRes = doAssign(dRes, sRes, hCfg);
        nClusters = numel(sRes.clusterCenters);
    end

    % remove small clusters
    for iRepeat = 1:1000
        sRes = removeSmall(dRes, sRes, hCfg);
        if nClusters == numel(sRes.clusterCenters)
            break;
        end
        sRes = doAssign(dRes, sRes, hCfg);
        nClusters = numel(sRes.clusterCenters);
    end

    hCfg.updateLog('assignClusters', sprintf('Finished initial assignment (%d clusters)', nClusters), 0, 1);
end

%% LOCAL FUNCTIONS
function sRes = computeCenters(dRes, sRes, hCfg)
    %COMPUTECENTERS Find cluster centers
    if ~isfield(dRes, 'spikesBySite')
        dRes.spikesBySite = arrayfun(@(iSite) dRes.spikes(dRes.spikeSites==iSite), hCfg.siteMap, 'UniformOutput', 0);
    end

    if strcmp(hCfg.RDDetrendMode, 'local')      % perform detrending site by site
        sRes.clusterCenters = jrclust.sort.detrendRhoDelta(sRes, dRes.spikesBySite, 1, hCfg);
    elseif strcmp(hCfg.RDDetrendMode, 'global') % detrend over all sites
        sRes.clusterCenters = jrclust.sort.detrendRhoDelta(sRes, dRes.spikesBySite, 0, hCfg);
    elseif strcmp(hCfg.RDDetrendMode, 'logz')   % identify centers with sufficiently high z-scores
        % sRes.clusterCenters = log_ztran_(sRes.spikeRho, sRes.spikeDelta, hCfg.log10RhoCut, 4 + hCfg.log10DeltaCut);
        x = log10(sRes.spikeRho(:));
        y = log10(sRes.spikeDelta(:));

        mask = isfinite(x) & isfinite(y);
        y(mask) = jrclust.utils.zscore(y(mask));

        sRes.clusterCenters = find(x >= hCfg.log10RhoCut & y >= 4 + hCfg.log10DeltaCut);
    elseif strcmp(hCfg.RDDetrendMode, 'regress') % regression method
        sRes.clusterCenters = jrclust.sort.regressCenters(sRes, dRes.spikesBySite, hCfg.log10DeltaCut);
    else                                            % don't detrend
        sRes.clusterCenters = find(sRes.spikeRho(:) > 10^(hCfg.log10RhoCut) & sRes.spikeDelta(:) > 10^(hCfg.log10DeltaCut));
    end
end

function sRes = doAssign(dRes, sRes, hCfg)
    %DOASSIGN Assign spikes to clusters
    nSpikes = numel(sRes.ordRho);
    nClusters = numel(sRes.clusterCenters);

    if isempty(sRes.spikeClusters)
        sRes.spikeClusters = zeros([nSpikes, 1], 'int32');
        sRes.spikeClusters(sRes.clusterCenters) = 1:nClusters;
    end

    % one or no center, assign all spikes to one cluster
    if numel(sRes.clusterCenters) == 0 || numel(sRes.clusterCenters) == 1
        sRes.spikeClusters = ones([nSpikes, 1], 'int32');
        sRes.clusterCenters = sRes.ordRho(1);
    else
        unassigned = sRes.spikeClusters <= 0;
        canAssign =  sRes.spikeClusters(sRes.spikeNeigh) > 0;
        doAssign = unassigned & canAssign;

        while any(doAssign)
            hCfg.updateLog('assignIter', sprintf('%d/%d spikes unassigned, %d can be assigned', ...
                                                 sum(unassigned), nSpikes, sum(doAssign)), 0, 0);
            sRes.spikeClusters(doAssign) = sRes.spikeClusters(sRes.spikeNeigh(doAssign));

            unassigned = sRes.spikeClusters <= 0;
            canAssign  = sRes.spikeClusters(sRes.spikeNeigh) > 0;
            doAssign = unassigned & canAssign;
        end
    end

    nClusters = numel(sRes.clusterCenters);

    % count spikes in clusters
    sRes.spikesByCluster = arrayfun(@(iC) find(sRes.spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
    sRes.unitCount = cellfun(@numel, sRes.spikesByCluster);
    sRes.clusterSites = double(arrayfun(@(iC) mode(dRes.spikeSites(sRes.spikesByCluster{iC})), 1:nClusters));
end

function sRes = splitDist(dRes, sRes, hCfg)
    nClusters = numel(sRes.clusterCenters);
    for jCluster = 1:nClusters
        % get the number of unique sites for this cluster
        jSpikes = sRes.spikesByCluster{jCluster};
        jSites = dRes.spikeSites(jSpikes);
        uniqueSites = unique(jSites);
        if numel(unique(jSites)) == 1
            continue;
        end

        % order spike sites by density, descending
        jRho = sRes.spikeRho(jSpikes);
        [~, ordering] = sort(jRho, 'descend');
        jSpikes = jSpikes(ordering);
        jSites = jSites(ordering);

        siteOrdering = arrayfun(@(k) find(jSites == k, 1), uniqueSites);
        [~, siteOrdering] = sort(siteOrdering);
        uniqueSites = uniqueSites(siteOrdering);

        siteLocs = hCfg.siteLoc(uniqueSites, :);
        siteDists = pdist2(siteLocs, siteLocs);

        % find pairwise distances which exceed the merge radius
        isFar = siteDists > 2*hCfg.evtMergeRad;
        [r, c] = find(isFar);
        if isempty(r)
            continue;
        end

        splitOff = jSpikes(arrayfun(@(k) find(jSites == k, 1), unique(uniqueSites(r(r > c)))));
        sRes.clusterCenters = [sRes.clusterCenters; splitOff];
    end
end

function sRes = removeSmall(dRes, sRes, hCfg)
    hCfg.minClusterSize = max(hCfg.minClusterSize, 2*size(dRes.spikeFeatures, 1));

    % remove small clusters
    smallClusters = sRes.unitCount <= hCfg.minClusterSize;

    if any(smallClusters)
        sRes.clusterCenters(smallClusters) = [];
        sRes.spikeClusters = [];
    end
end
