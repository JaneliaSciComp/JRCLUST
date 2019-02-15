% function positions = spk_pos_(S0, trFet_spk)
function positions = spikePos(spikeSites, spikeFeatures, hCfg)
    %SPIKEPOS
    nSites_spk = 1 + hCfg.nSiteDir*2 - hCfg.nSitesExcl;

    mrVp = squeeze(spikeFeatures(1:nSites_spk, 1, :)) .^ 2;
    vrVp = sum(mrVp);

    miSites_spk = single(hCfg.siteNeighbors(1:nSites_spk, spikeSites));
    mrX_spk = reshape(hCfg.siteLoc(miSites_spk,1), size(miSites_spk));
    mrY_spk = reshape(hCfg.siteLoc(miSites_spk,2), size(miSites_spk));

    positions = zeros(numel(spikeSites), 2, 'single');
    positions(:, 1) = sum(mrVp .* mrX_spk) ./ vrVp;
    positions(:, 2) = sum(mrVp .* mrY_spk) ./ vrVp;
end
