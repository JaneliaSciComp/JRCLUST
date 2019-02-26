function rmOutlierSpikes(obj)
    %RMOUTLIERSPIKES Mahalanobis-distance based outlier removal
    if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

    if obj.hCfg.outlierThresh == 0
        return;
    end

    % This is troubling:
    % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
    % > In mahal (line 49)
    % TODO: investigate
    warning off;
    for iCluster = 1:obj.nClusters
        iSite = obj.clusterSites(iCluster);
        iSpikes = obj.spikesByCluster{iCluster};

        if isempty(iSpikes)
            continue;
        end

        iFeatures = squeeze(obj.spikeFeatures(:, 1, iSpikes));

        % if iSite is secondary site for any spikes in this cluster, prefer
        % features computed from those over primary spike features
        if size(obj.spikeFeatures, 2) >= 2
            iFeatures2 = squeeze(obj.spikeFeatures(:, 2, iSpikes));
            onSite2 = find(obj.spikeSites2(iSpikes) == iSite);
            iFeatures(:, onSite2) = iFeatures2(:, onSite2);
        end

        iFeatures = iFeatures';

        try % MAD transform of log self-Mahalanobis distance
            lastwarn(''); % reset last warning to catch it

            iDist = jrclust.utils.madScore(log(mahal(iFeatures, iFeatures)));
            [wstr, wid] = lastwarn();
            if strcmp(wid, 'MATLAB:nearlySingularMatrix')
                error(wstr);
            end
        catch
            continue;
        end

        exclude = (iDist > obj.hCfg.outlierThresh);
        if any(exclude)
            obj.spikesByCluster{iCluster} = iSpikes(~exclude);
            obj.spikeClusters(iSpikes(exclude)) = 0; % classify as noise
            obj.unitCount(iCluster) = numel(obj.spikesByCluster{iCluster});
        end
    end
    warning on;
end