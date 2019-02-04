function sRes = assignClusters(dRes, sRes, hCfg)
    %ASSIGNCLUSTERS Given rho-delta information, assign spikes to clusters
    sRes = computeCenters(dRes, sRes, hCfg);
    sRes.spikeClusters = [];
    sRes = doAssignClusters(dRes, sRes, hCfg);
end

%% LOCAL FUNCTIONS
function sRes = computeCenters(dRes, sRes, hCfg)
    %COMPUTECENTERS Find cluster centers
    if ~isfield(dRes, 'spikesBySite')
        dRes.spikesBySite = arrayfun(@(iSite) dRes.spikes(dRes.spikeSites==iSite), hCfg.siteMap, 'UniformOutput', 0);
    end

    if strcmp(hCfg.RDDetrendMode, 'local')      % perform detrending site by site
        sRes.clusterCenters = jrclust.cluster.densitypeaks.detrendRhoDelta(sRes, dRes.spikesBySite, 1, hCfg);
    elseif strcmp(hCfg.RDDetrendMode, 'global') % detrend over all sites
        sRes.clusterCenters = jrclust.cluster.densitypeaks.detrendRhoDelta(sRes, dRes.spikesBySite, 0, hCfg);
    elseif strcmp(hCfg.RDDetrendMode, 'logz')   % identify centers with sufficiently high z-scores
        % sRes.clusterCenters = log_ztran_(sRes.spikeRho, sRes.spikeDelta, hCfg.log10RhoCut, 4 + hCfg.log10DeltaCut);
        x = log10(sRes.spikeRho(:));
        y = log10(sRes.spikeDelta(:));

        mask = isfinite(x) & isfinite(y);
        y(mask) = jrclust.utils.zscore(y(mask));

        sRes.clusterCenters = find(x >= hCfg.log10RhoCut & y >= 4 + hCfg.log10DeltaCut);
    elseif strcmp(hCfg.RDDetrendMode, 'regress') % regression method
        sRes.clusterCenters = jrclust.cluster.densitypeaks.regressCenters(sRes, dRes.spikesBySite, hCfg.log10DeltaCut);
    else                                            % don't detrend
        sRes.clusterCenters = find(sRes.spikeRho(:) > 10^(hCfg.log10RhoCut) & sRes.spikeDelta(:) > 10^(hCfg.log10DeltaCut));
    end
end

function sRes = doAssignClusters(dRes, sRes, hCfg)
    nRepeatMax = 1000;
    if isempty(sRes.spikeClusters)
        nClustersPrev = [];
    else
        nClustersPrev = sRes.nClusters;
    end

    removedClusters = 0;
    if hCfg.verbose
        fprintf('assigning clusters, nClusters:%d\n', numel(sRes.clusterCenters));
        t = tic;
    end

    % assign spikes to clusters
    for iRepeat = 1:nRepeatMax % repeat 1000 times max
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
            nNeigh = sRes.spikeNeigh(sRes.ordRho);

            for i = 1:10
                unassigned = find(sRes.spikeClusters(sRes.ordRho) <= 0);
                if isempty(unassigned)
                    break;
                end

                unassigned = unassigned(:)';

                for j = unassigned
                    sRes.spikeClusters(sRes.ordRho(j)) = sRes.spikeClusters(nNeigh(j));
                end
                nUnassigned = sum(sRes.spikeClusters <= 0);

                if nUnassigned == 0
                    break;
                end

                if hCfg.verbose
                    fprintf('i: %d, n0 = %d, ', i, nUnassigned);
                end
            end
            sRes.spikeClusters(sRes.spikeClusters <= 0) = 1; %background
        end

        hCfg.minClusterSize = max(hCfg.minClusterSize, 2*size(dRes.spikeFeatures, 1));

        % count spikes in clusters
        sRes.spikesByCluster = arrayfun(@(iC) find(sRes.spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
        sRes.unitCount = cellfun(@numel, sRes.spikesByCluster);
        sRes.clusterSites = double(arrayfun(@(iC) mode(dRes.spikeSites(sRes.spikesByCluster{iC})), 1:nClusters));

        % remove small clusters
        smallClusters = find(sRes.unitCount <= hCfg.minClusterSize);
        if isempty(smallClusters) % done!
            break;
        end

        % still here? try again
        sRes.clusterCenters(smallClusters) = [];
        sRes.spikeClusters = [];
        removedClusters = removedClusters + numel(smallClusters);

        if iRepeat == nRepeatMax
            warning('assignClusters: exceeded nRepeatMax = %d\n', nRepeatMax);
        end
    end % for

    if hCfg.verbose
        fprintf('\n\ttook %0.1fs. Removed %d clusters having <%d spikes: %d->%d\n', toc(t), removedClusters, hCfg.minClusterSize, nClustersPrev, nClusters);
    end
end
