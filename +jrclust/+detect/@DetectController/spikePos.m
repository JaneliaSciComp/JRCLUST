function positions = spikePos(obj, spikeSites, spikeFeatures)
    %SPIKEPOS
    nSitesEvt = 1 + obj.hCfg.nSiteDir*2 - obj.hCfg.nSitesExcl;

    featuresSquared = squeeze(spikeFeatures(1:nSitesEvt, 1, :)) .^ 2;
    vrVp = sum(featuresSquared);

    spikeNeighbors = single(obj.hCfg.siteNeighbors(1:nSitesEvt, spikeSites));
    spikeX = reshape(obj.hCfg.siteLoc(spikeNeighbors, 1), size(spikeNeighbors));
    spikeY = reshape(obj.hCfg.siteLoc(spikeNeighbors, 2), size(spikeNeighbors));

    positions = zeros(numel(spikeSites), 2, 'single');
    positions(:, 1) = sum(featuresSquared .* spikeX) ./ vrVp;
    positions(:, 2) = sum(featuresSquared .* spikeY) ./ vrVp;
end
