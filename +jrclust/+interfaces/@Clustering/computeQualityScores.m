function computeQualityScores(obj, updateMe)
    %COMPUTEQUALITYSCORES Get cluster quality scores
    if nargin < 2
        updateMe = [];
    end

    obj.hCfg.updateLog('qualScores', 'Computing cluster quality scores', 1, 0);

    unitVmin = squeeze(min(obj.meanWfGlobal));
    unitVmax = squeeze(max(obj.meanWfGlobal));
    unitVminRaw = squeeze(min(obj.meanWfGlobalRaw));
    unitVmaxRaw = squeeze(max(obj.meanWfGlobalRaw));

    unitVpp_ = jrclust.utils.rowColSelect(unitVmax - unitVmin, obj.clusterSites, 1:obj.nClusters);
    unitVppRaw_ = jrclust.utils.rowColSelect(unitVmaxRaw - unitVminRaw, obj.clusterSites, 1:obj.nClusters);
    unitPeaksRaw_ = jrclust.utils.rowColSelect(unitVminRaw, obj.clusterSites, 1:obj.nClusters);

    % compute unitIsoDist_, unitLRatio_, unitISIRatio_
    nSamples2ms = round(obj.hCfg.sampleRate * .002);
    nSamples20ms = round(obj.hCfg.sampleRate * .02);

    if isempty(updateMe)
        [unitIsoDist_, unitLRatio_, unitISIRatio_] = deal(nan(obj.nClusters, 1));
        updateMe = 1:obj.nClusters;
    else
        unitIsoDist_ = obj.unitIsoDist;
        unitLRatio_ = obj.unitLRatio;
        unitISIRatio_ = obj.unitISIRatio;
        updateMe = updateMe(:)';
    end

    % This is troubling:
    % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
    % > In mahal (line 49)
    % TODO: investigate
    warning off;
    for iCluster = updateMe
        clusterSpikes_ = obj.spikesByCluster{iCluster};
        % Compute ISI ratio
        clusterTimes_ = obj.spikeTimes(clusterSpikes_);
        diffCtimes = diff(clusterTimes_);
        
        % define ISI ratio as #(ISI <= 2ms)/#(ISI <= 20ms)
        unitISIRatio_(iCluster) = sum(diffCtimes <= nSamples2ms)./sum(diffCtimes <= nSamples20ms);

        % Compute L-ratio and isolation distance (use neighboring features)
        iSite = obj.clusterSites(iCluster);

        % find spikes whose primary or secondary spikes live on iSite
        iSite1Spikes = obj.spikesBySite{iSite};
        if isprop(obj, 'spikesBySite2') && ~isempty(obj.spikesBySite2{iSite})
            iSite2Spikes = obj.spikesBySite2{iSite};
            siteFeatures = [squeeze(obj.spikeFeatures(:, 1, iSite1Spikes)), squeeze(obj.spikeFeatures(:, 2, iSite2Spikes))];
            iSiteSpikes = [iSite1Spikes(:); iSite2Spikes(:)];
        else
            siteFeatures = squeeze(obj.spikeFeatures(:, 1, iSite1Spikes));
            iSiteSpikes = iSite1Spikes(:);
        end

        isOnSite = (obj.spikeClusters(iSiteSpikes) == iCluster);
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

    obj.unitISIRatio = unitISIRatio_;
    obj.unitIsoDist = unitIsoDist_;
    obj.unitLRatio = unitLRatio_;
    obj.unitPeaksRaw = unitPeaksRaw_; % unitPeaks is set elsewhere
    obj.unitVpp = unitVpp_;
    obj.unitVppRaw = unitVppRaw_;

    obj.hCfg.updateLog('qualScores', 'Finished computing cluster quality scores', 0, 1);
end