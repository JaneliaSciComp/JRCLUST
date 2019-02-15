function clusterMean = nanmeanInt16(spikeWindows, iSite, sites, hCfg)
    %NANMEANINT16 
    iSiteNeighbors = hCfg.siteNeighbors(:, iSite);
    trWav = nan([size(spikeWindows, 1), numel(iSiteNeighbors), numel(sites)], 'single');

    uniqueSites = unique(sites);
    nUniqueSites = numel(uniqueSites);
    uniqueNeighbors = hCfg.siteNeighbors(:, uniqueSites);

    for jSite = 1:nUniqueSites
        iSiteUnique = uniqueSites(jSite);
        viSpk_ = find(sites == iSiteUnique);

        [~, viSite1a_, viSite1b_] = intersect(iSiteNeighbors, uniqueNeighbors(:, jSite));
        if isempty(viSite1a_)
            continue;
        end

        trWav(:, viSite1a_, viSpk_) = spikeWindows(:, viSite1b_, viSpk_);
    end

    clusterMean = nanmean(trWav, 3);
    clusterMean = jrclust.utils.meanSubtract(clusterMean); %122717 JJJ
end