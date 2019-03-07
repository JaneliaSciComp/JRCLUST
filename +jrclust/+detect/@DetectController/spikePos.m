function positions = spikePos(obj, spikeSites, spikeFeatures)
    %SPIKEPOS Compute a feature-weighted position for each spike
    nSitesEvt = 1 + obj.hCfg.nSiteDir*2 - obj.hCfg.nSitesExcl;

    featuresSquared = squeeze(spikeFeatures(1:nSitesEvt, 1, :)) .^ 2;
    if iscolumn(featuresSquared)
        featuresSquared = featuresSquared';
    end

    vrVp = sum(featuresSquared, 1);

    spikeNeighbors = single(obj.hCfg.siteNeighbors(1:nSitesEvt, spikeSites));
    spikeX = reshape(obj.hCfg.siteLoc(spikeNeighbors, 1), size(spikeNeighbors));
    spikeY = reshape(obj.hCfg.siteLoc(spikeNeighbors, 2), size(spikeNeighbors));

    positions = zeros(numel(spikeSites), 2, 'single');
    positions(:, 1) = sum(featuresSquared .* spikeX) ./ vrVp;
    positions(:, 2) = sum(featuresSquared .* spikeY) ./ vrVp;
end
