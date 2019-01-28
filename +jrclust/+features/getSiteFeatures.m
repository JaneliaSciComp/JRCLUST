% function [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk)
function [siteFeatures, spikes, n1, n2, spikeOrder] = getSiteFeatures(spikeFeatures, site, spikeData, hCfg)
    %GETSITEFEATURES Get features occurring on primary and secondary (optionally tertiary) sites
    [siteFeatures, spikes, n1, n2, spikeOrder] = deal([]);

    nPeaksFeatures = hCfg.nPeaksFeatures; % or maybe size(spikeFeatures, 2) ? 
    timeFeatureFactor = hCfg.getOr('timeFeatureFactor', 0); % TW

    if isfield(spikeData, 'spikes1') && ~isempty(spikeData.spikes1)
        spikes1 = int32(spikeData.spikes1);
    else
        return
    end

    if isfield(spikeData, 'spikes2') && ~isempty(spikeData.spikes2)
        spikes2 = int32(spikeData.spikes2);
    else
        spikes2 = [];
        nPeaksFeatures = 1;
    end

    if isfield(spikeData, 'spikes3') && ~isempty(spikeData.spikes3)
        spikes3 = int32(spikeData.spikes3);
    else
        spikes3 = [];
    end

    if isfield(spikeData, 'vlRedo_spk') && ~isempty(spikeData.vlRedo_spk)
        spikes1 = spikes1(spikeData.vlRedo_spk(spikes1));
        spikes2 = spikes2(spikeData.vlRedo_spk(spikes2));
        spikes3 = spikes3(spikeData.vlRedo_spk(spikes3));
    end

    n1 = numel(spikes1);
    n2 = numel(spikes2);

    % get features for each site
    if nPeaksFeatures == 1
        siteFeatures = [squeeze(spikeFeatures(:, 1, spikes1)); single(spikeData.spikeTimes(spikes1))']; % TW
        siteFeatures(end, :) = timeFeatureFactor*std(siteFeatures(1, :)).*siteFeatures(end, :)./std(siteFeatures(end, :)); % TW

        spikes = spikes1;
    elseif nPeaksFeatures == 2
        siteFeatures = [squeeze(spikeFeatures(:, 1, spikes1)), squeeze(spikeFeatures(:, 2, spikes2)); single(spikeData.spikeTimes([spikes1; spikes2]))']; % TW
        siteFeatures(end, :) = timeFeatureFactor*std(siteFeatures(1, :)).*siteFeatures(end, :)./std(siteFeatures(end, :)); % TW

        spikes = [spikes1; spikes2];
    else % nPeaksFeatures == 3
        siteFeatures = [squeeze(spikeFeatures(:, 1, spikes1)), squeeze(spikeFeatures(:, 2, spikes2)), squeeze(spikeFeatures(:, 3, spikes3)); single(spikeData.spikeTimes([spikes1; spikes2; spikes3]))'];
        siteFeatures(end, :) = timeFeatureFactor*std(siteFeatures(1, :)).*siteFeatures(end, :)./std(siteFeatures(end, :)); % TW

        spikes = [spikes1; spikes2; spikes3];
        n2 = n2 + numel(spikes3);
    end

    % weight features by distance from site, greater for nearer sites
    try
        nSites = hCfg.nSitesEvt;

        if hCfg.weightFeatures && nSites >= hCfg.minSitesWeightFeatures
            nFeaturesPerSite = size(siteFeatures, 1) / nSites;

            weights = distWeight(site, nSites, hCfg);
            weights = repmat(weights(:), [nFeaturesPerSite, 1]);

            siteFeatures = bsxfun(@times, siteFeatures, weights(:));
        end
    catch ME
        warning('error in distWeight: ''%s'', using unweighted features', ME.message);
    end

    spikeOrder = jrclust.utils.rankorder(spikes, 'ascend');
end

%% LOCAL FUNCTIONS
function weights = distWeight(site, nSites, hCfg)
    %DISTWEIGHT Compute weights for neighboring sites, larger for nearer
    siteLoc = hCfg.siteLoc(hCfg.siteNeighbors(1:nSites, site), :);
    siteDists = pdist2(siteLoc(1, :), siteLoc);

    weights = 2.^(-siteDists/hCfg.evtDetectRad);
    weights = weights(:);
end
