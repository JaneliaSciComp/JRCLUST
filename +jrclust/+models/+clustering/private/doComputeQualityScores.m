function scores = doComputeQualityScores(hClust, updateMe)
    %DOCOMPUTEQUALSCORES Compute cluster quality scores
    if nargin < 2
        updateMe = [];
    end

    if hClust.hCfg.verbose
        t1 = tic;
        fprintf('Calculating cluster quality...\n');
    end

    unitVmin = squeeze(min(hClust.meanWfGlobal));
    unitVmax = squeeze(max(hClust.meanWfGlobal));
    unitVminRaw = squeeze(min(hClust.meanWfGlobalRaw));
    unitVmaxRaw = squeeze(max(hClust.meanWfGlobalRaw));

    unitVpp_ = jrclust.utils.rowColSelect(unitVmax - unitVmin, hClust.clusterSites, 1:hClust.nClusters);
    unitPeaks_ = jrclust.utils.rowColSelect(unitVmin, hClust.clusterSites, 1:hClust.nClusters);
    unitVppRaw_ = jrclust.utils.rowColSelect(unitVmaxRaw - unitVminRaw, hClust.clusterSites, 1:hClust.nClusters);
    unitPeaksRaw_ = jrclust.utils.rowColSelect(unitVminRaw, hClust.clusterSites, 1:hClust.nClusters);

    try
        siteRMS_ = jrclust.utils.bit2uV(single(hClust.siteThresh(:))/hClust.hCfg.qqFactor, hClust.hCfg);
        unitSNR_ = abs(unitPeaks_)./siteRMS_(hClust.clusterSites);
        nSitesOverThresh_ = sum(bsxfun(@lt, unitVmin, - siteRMS_*hClust.hCfg.qqFactor), 1)';
    catch ME
        [siteRMS_, unitSNR_, nSitesOverThresh_] = deal([]);
        warning('RMS, SNR, nSitesOverThresh not set: %s', ME.message);
    end

    % compute unitIsoDist_, unitLRatio_, unitISIRatio_
    nSamples2ms = round(hClust.hCfg.sampleRate * .002);
    nSamples20ms = round(hClust.hCfg.sampleRate * .02);

    if isempty(updateMe)
        [unitIsoDist_, unitLRatio_, unitISIRatio_] = deal(nan(hClust.nClusters, 1));
        updateMe = 1:hClust.nClusters;
    else
        unitIsoDist_ = hClust.unitIsoDist;
        unitLRatio_ = hClust.unitLRatio;
        unitISIRatio_ = hClust.unitISIRatio;
        updateMe = updateMe(:)';
    end

    % This is troubling:
    % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
    % > In mahal (line 49)
    % TODO: investigate
    warning off;
    for iCluster = updateMe
        clusterSpikes_ = hClust.spikesByCluster{iCluster};
        % Compute ISI ratio
        clusterTimes_ = hClust.spikeTimes(clusterSpikes_);
        diffCtimes = diff(clusterTimes_);
        
        % define ISI ratio as #(ISI <= 2ms)/#(ISI <= 20ms)
        unitISIRatio_(iCluster) = sum(diffCtimes <= nSamples2ms)./sum(diffCtimes <= nSamples20ms);

        % Compute L-ratio and isolation distance (use neighboring features)
        iSite = hClust.clusterSites(iCluster);

        % find spikes whose primary or secondary spikes live on iSite
        iSite1Spikes = hClust.spikesBySite{iSite};
        iSite2Spikes = hClust.spikesBySite2{iSite};

        if isempty(iSite2Spikes)
            siteFeatures = squeeze(hClust.spikeFeatures(:, 1, iSite1Spikes));
            iSiteSpikes = iSite1Spikes(:);
        else
            siteFeatures = [squeeze(hClust.spikeFeatures(:, 1, iSite1Spikes)), squeeze(hClust.spikeFeatures(:, 2, iSite2Spikes))];
            iSiteSpikes = [iSite1Spikes(:); iSite2Spikes(:)];
        end

        isOnSite = (hClust.spikeClusters(iSiteSpikes) == iCluster);
        nSpikesOnSite = sum(isOnSite);

        [unitLRatio_(iCluster), unitIsoDist_(iCluster)] = deal(nan);

        try
            lastwarn(''); % reset last warning to catch it
            mDists = mahal(siteFeatures', siteFeatures(:, isOnSite)');

            [wstr, wid] = lastwarn();
            if strcmp(wid, 'MATLAB:nearlySingularMatrix')
                error(wstr);
            end
        catch
            continue;
        end

        mDistsOutside = mDists(~isOnSite);
        if isempty(mDistsOutside)
            continue;
        end

        unitLRatio_(iCluster) = sum(1 - chi2cdf(mDistsOutside, size(siteFeatures,2)))/nSpikesOnSite;

        % compute isolation distance
        if mean(isOnSite) > .5
            continue;
        end

        sorted12 = sort(mDistsOutside);
        unitIsoDist_(iCluster) = sorted12(nSpikesOnSite);
    end % for
    warning on;

    scores = struct();
    scores.nSitesOverThresh = nSitesOverThresh_;
    scores.siteRMS = siteRMS_;
    scores.unitISIRatio = unitISIRatio_;
    scores.unitIsoDist = unitIsoDist_;
    scores.unitLRatio = unitLRatio_;
    scores.unitPeaksRaw = unitPeaksRaw_; % unitPeaks is set elsewhere
    scores.unitSNR = unitSNR_;
    scores.unitVpp = unitVpp_;
    scores.unitVppRaw = unitVppRaw_;

    if hClust.hCfg.verbose
        fprintf('\ttook %0.1fs\n', toc(t1));
    end
end
