function [hClust, S0] = autoMerge(hClust, hCfg, doAssign)
    %AUTOMERGE Automatically merge clusters
    if nargin < 4
        doAssign = true;
    end

%     if doAssign
%         S_clu = postCluster_(S_clu, P);
%     end

%     hCluster = rmfield_(hCluster, 'tmrWav_clu', 'vrPosX_clu', 'vrPosY_clu', 'vrVpp_clu', ...
%         'vrVpp_uv_clu', 'vrVmin_uv_clu', 'vrSnr_clu', 'vnSite_clu', ...
%         'vrIsoDist_clu', 'vrLRatio_clu', 'vrIsiRatio_clu');
    hClust.refresh();
    hClust.orderClusters('clusterSites');
    hClust.clearNotes();

    if hCfg.outlierThresh > 0
        hClust = rmOutlierSpikes(hClust, hCfg.outlierThresh);
    end

    hClust = post_merge_wav_(hClust, spikeData, hCfg);

    hClust.refresh();
    hClust.orderClusters('clusterSites');
    hClust = S_clu_update_wav_(hClust, hCfg);

    % set diagonal element
    [hClust, S0] = S_clu_commit_(hClust, 'jrclust.clustering.autoMerge');
    hClust.mrWavCor = set_diag_(hClust.mrWavCor, S_clu_self_corr_(hClust, [], S0));
    hClust.P = hCfg;
    hClust = S_clu_position_(hClust);
    hClust.csNote_clu = cell(hClust.nClu, 1); %reset note
    hClust = S_clu_quality_(hClust, hCfg);
    [hClust, S0] = S_clu_commit_(hClust, 'jrclust.clustering.autoMerge');
end

%% LOCAL FUNCTIONS
function hClust = rmOutlierSpikes(hClust, thresh)
    %RMOUTLIERSPIKES Mahalanobis-distance based outlier removal
    for iCluster = 1:hClust.nClusters
        iSite = hClust.clusterSites(iCluster);
        iSpikes = hClust.clusterSpikes{iCluster};

        if isempty(iSpikes)
            continue;
        end

        iFeatures = squeeze(hClust.spikeFeatures(:, 1, iSpikes));

        % if iSite is secondary site for any spikes in this cluster, prefer
        % features computed from those over primary spike features
        if size(hClust.spikeFeatures, 2) >= 2
            iFeatures2 = squeeze(hClust.spikeFeatures(:, 2, iSpikes));
            onSite2 = find(hClust.siteSpikes2(iSpikes) == iSite);
            iFeatures(:, onSite2) = iFeatures2(:, onSite2);
        end

        iFeatures = iFeatures';

        try % MAD transform of log self-Mahalanobis distance
            iDist = jrclust.utils.madScore(log(mahal(iFeatures, iFeatures)));
        catch
            continue;
        end

        exclude = (iDist > thresh);
        if any(exclude)
            hClust.clusterSpikes{iCluster} = iSpikes(~exclude);
            hClust.spikeClusters(iSpikes(exclude)) = 0; % classify as noise
            hClust.clusterCounts(iCluster) = numel(hClust.clusterSpikes{iCluster});
        end

        fprintf('.');
    end
end
