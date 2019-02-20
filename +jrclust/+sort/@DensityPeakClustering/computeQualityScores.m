function computeQualityScores(obj, updateMe)
    %COMPUTEQUALITYSCORES Get cluster quality scores
    computeQualityScores@jrclust.interfaces.Clustering(obj, updateMe);

    unitVmin = squeeze(min(obj.meanWfGlobal));
    unitPeaks_ = jrclust.utils.rowColSelect(unitVmin, obj.clusterSites, 1:obj.nClusters);

    try
        siteRMS_ = jrclust.utils.bit2uV(single(obj.siteThresh(:))/obj.hCfg.qqFactor, obj.hCfg);
        unitSNR_ = abs(unitPeaks_)./siteRMS_(obj.clusterSites);
        nSitesOverThresh_ = sum(bsxfun(@lt, unitVmin, - siteRMS_*obj.hCfg.qqFactor), 1)';
    catch ME
        [siteRMS_, unitSNR_, nSitesOverThresh_] = deal([]);
        warning('RMS, SNR, nSitesOverThresh not set: %s', ME.message);
    end

    obj.unitSNR = unitSNR_;
    obj.nSitesOverThresh = nSitesOverThresh_;
    obj.siteRMS = siteRMS_;
end

