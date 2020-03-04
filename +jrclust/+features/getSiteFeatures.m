% function [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk)
function [siteFeatures, spikes, n1, n2, spikeOrder] = getSiteFeatures(spikeFeatures, site, spikeData, hCfg)
    %GETSITEFEATURES Get features occurring on primary and secondary (optionally tertiary) sites
    [siteFeatures, spikes, n1, n2, spikeOrder] = deal([]);
    [spikes2, spikes3, timeFeature] = deal([]);

    [nFeatures, nPeaksFeatures, ~] = size(spikeFeatures);
    timeFeatureFactor = hCfg.getOr('timeFeatureFactor', 0); % TW

    if ~isfield(spikeData, 'spikes1') || isempty(spikeData.spikes1)
        return;
    end

    spikes1 = int32(spikeData.spikes1);

    if isfield(spikeData, 'spikes2')
        spikes2 = int32(spikeData.spikes2);
    end

    if isfield(spikeData, 'spikes3')
        spikes3 = int32(spikeData.spikes3);
    end

    if isfield(spikeData, 'vlRedo_spk') && ~isempty(spikeData.vlRedo_spk)
        spikes1 = spikes1(spikeData.vlRedo_spk(spikes1));
        spikes2 = spikes2(spikeData.vlRedo_spk(spikes2));
        spikes3 = spikes3(spikeData.vlRedo_spk(spikes3));
    end

    % compute time feature (last entry of feature vector
    if timeFeatureFactor ~= 0
        times1 = single(spikeData.spikeTimes(spikes1))'; % row vector
        times2 = single(spikeData.spikeTimes(spikes2))'; % row vector or empty
        times3 = single(spikeData.spikeTimes(spikes3))'; % row vector or empty
        timeFeature = timeFeatureFactor * [times1 times2 times3];
    end

    % features from the first peak
    n1 = numel(spikes1);
    siteFeatures = getPeakFeature(spikeFeatures(:, 1, spikes1), nFeatures, n1);

    % features from the second peak
    if nPeaksFeatures > 1
        n2 = numel(spikes2);
        sf2 = getPeakFeature(spikeFeatures(:, 2, spikes2), nFeatures, n2);
        siteFeatures = [siteFeatures sf2];
    end

    % features from the third peak
    if nPeaksFeatures > 2
        n3 = numel(spikes3);
        sf3 = getPeakFeature(spikeFeatures(:, 3, spikes3), nFeatures, n3);
        siteFeatures = [siteFeatures sf3];
        n2 = n2 + n3;
    end

    % scale time feature (no effect if empty) and add to siteFeatures
    timeFeature = timeFeature*std(siteFeatures(1, :))/std(timeFeature);
    siteFeatures = [siteFeatures; timeFeature];

    % weight features by distance from site, greater for nearer sites
    nSites = hCfg.nSitesEvt;

    try
        if hCfg.weightFeatures && nSites >= hCfg.minSitesWeightFeatures
            nFeaturesPerSite = size(siteFeatures, 1) / nSites;

            weights = distWeight(site, nSites, hCfg);
            weights = repmat(weights(:), [nFeaturesPerSite, 1]);

            siteFeatures = bsxfun(@times, siteFeatures, weights(:));
        end
    catch ME
        warning('error in distWeight: ''%s'', using unweighted features', ME.message);
    end

    % concatenate all spikes
    spikes = [spikes1; spikes2; spikes3];
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

function feat = getPeakFeature(feat, nFeatures, nSpikes)
    %GETPEAKFEATURES Reshape feat to nFeatures x nSpikes if necessary.
    feat = squeeze(feat);
    [d1, d2] = size(feat);

    % either 1, in which case a scalar, or n > 1, in which case leave alone
    if nSpikes == nFeatures || (d1 == nFeatures && d2 == nSpikes)
        return;
    end

    if d1 == nSpikes && d2 == nFeatures
        feat = feat';
    else
        error('Expected dimensions [%d, %d], got [%d, %d]', feat, nSpikes, d1, d2);
    end
end
