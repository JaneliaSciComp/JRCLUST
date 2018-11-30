%--------------------------------------------------------------------------
function hClust = S_clu_quality_(hClust, updateMe)
    % 7/5/17 JJJ: Added isolation distance, L-ratio,
    % TODO: update when deleting, splitting, or merging
    if nargin < 2
        updateMe = [];
    end

    t1 = tic;
    fprintf('Calculating cluster quality...\n');

    unitVmin = squeeze(min(hClust.meanWfGlobal));
    unitVmax = squeeze(max(hClust.meanWfGlobal));
    mrVmin_uv_clu = squeeze(min(hClust.meanWfGlobalRaw));
    mrVmax_uv_clu = squeeze(max(hClust.meanWfGlobalRaw));

    unitVpp_ = jrclust.utils.rowColSelect(unitVmax - unitVmin, hClust.clusterSites, 1:hClust.nClusters);
    unitPeaks_ = jrclust.utils.rowColSelect(unitVmin, hClust.clusterSites, 1:hClust.nClusters);
    unitVppRaw_ = jrclust.utils.rowColSelect(mrVmax_uv_clu-mrVmin_uv_clu, hClust.clusterSites, 1:hClust.nClusters);
    unitPeaksRaw_ = jrclust.utils.rowColSelect(mrVmin_uv_clu, hClust.clusterSites, 1:hClust.nClusters);

    try
        siteRMS_ = jrclust.utils.bit2uV(single(hClust.siteThresh(:))/hClust.hCfg.qqFactor, hClust.hCfg);
        unitSNR_ = abs(unitPeaks_) ./ siteRMS_(hClust.clusterSites);
        nSitesOverThresh_ = sum(bsxfun(@lt, unitVmin, - siteRMS_*hClust.hCfg.qqFactor), 1)';
    catch
        [siteRMS_, unitSNR_, nSitesOverThresh_] = deal([]);
        disp('no Sevt in memory.');
    end

    [unitIsoDist_, unitLRatio_, unitISIRatio_] = S_clu_quality2_(hClust, hClust.hCfg, updateMe);

    hClust.nSitesOverThresh = nSitesOverThresh_;
    hClust.siteRMS = siteRMS_;
    hClust.unitISIRatio = unitISIRatio_;
    hClust.unitIsoDist = unitIsoDist_;
    hClust.unitLRatio = unitLRatio_;
    hClust.unitPeaksRaw = unitPeaksRaw_; % unitPeaks set elsewhere
    hClust.unitSNR = unitSNR_;
    hClust.unitVpp = unitVpp_;
    hClust.unitVppRaw = unitVppRaw_;

    fprintf('\ttook %0.1fs\n', toc(t1));
end

function [unitIsoDist_, unitLRatio_, unitISIRatio_] = S_clu_quality2_(hClust, updateMe)
    % 7/5/17 JJJ: Berenyi2013 based
    nSamples_2ms = round(hClust.hCfg.sampleRate * .002);
    nSamples_20ms = round(hClust.hCfg.sampleRate * .02);

    if isempty(updateMe)
        [unitIsoDist_, unitLRatio_, unitISIRatio_] = deal(nan(hClust.nClusters, 1));
        updateMe = 1:hClust.nClusters;
    else
        unitIsoDist_ = hClust.unitIsoDist;
        unitLRatio_ = hClust.unitLRatio;
        unitISIRatio_ = hClust.unitISIRatio;
        updateMe = updateMe(:)';
    end

    for iClu = updateMe
        viSpk_clu1 = hClust.spikesBySite{iClu};
        % Compute ISI ratio
        viTime_clu1 = hClust.spikeTimes(viSpk_clu1);
        viDTime_clu1 = diff(viTime_clu1);
        unitISIRatio_(iClu) = sum(viDTime_clu1<=nSamples_2ms) ./ sum(viDTime_clu1<=nSamples_20ms);

        % Compute L-ratio an disodist (use neighboring features
        iSite_clu1 = hClust.viSite_clu(iClu);
        % find spikes whose primary or secondary spikes reside there
        viSpk1_local = hClust.spikesBySite{iSite_clu1}; %find(viSite_spk == iSite_clu1);
        viSpk2_local = hClust.spikesBySite2{iSite_clu1}; %find(viSite2_spk == iSite_clu1);
        if isempty(viSpk2_local)
            mrFet12_spk = squeeze_(hClust.spikeFeatures(:,1,viSpk1_local));
            viSpk12_local = viSpk1_local(:);
        else
            mrFet12_spk = [squeeze_(hClust.spikeFeatures(:,1,viSpk1_local)), squeeze_(hClust.spikeFeatures(:,2,viSpk2_local))]; %squeeze_(cat(3, hClust.spikeFeatures(:,1,viSpk1_local), hClust.spikeFeatures(:,2,viSpk2_local)))';
            viSpk12_local = [viSpk1_local(:); viSpk2_local(:)];
        end
        vlClu12_local = hClust.viClu(viSpk12_local) == iClu;
        nSpk_clu1 = sum(vlClu12_local);
        [unitLRatio_(iClu), unitIsoDist_(iClu)] = deal(nan);
        try
            vrMahal12 = mahal(mrFet12_spk', mrFet12_spk(:,vlClu12_local)');
        catch
            continue;
        end
        vrMahal12_out = vrMahal12(~vlClu12_local);
        if isempty(vrMahal12_out), continue; end
        unitLRatio_(iClu) = sum(1-chi2cdf(vrMahal12_out, size(mrFet12_spk,2))) / nSpk_clu1;

        % compute isolation distance
        if mean(vlClu12_local) > .5, continue; end
        sorted12 = sort(vrMahal12_out);
        unitIsoDist_(iClu) = sorted12(nSpk_clu1);
    end %for
end %func
