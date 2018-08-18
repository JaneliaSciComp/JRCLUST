%--------------------------------------------------------------------------
function mrPos_spk = spk_pos_(S0, spikeFeatures)
    if nargin<2, spikeFeatures = getSpikeFeatures(S0.P); end
    P = S0.P;
    nSites_spk = 1 + P.maxSite*2 - P.nSites_ref;
    %nSites_spk = size(spikeFeatures, 1);

    mrVp = squeeze_(spikeFeatures(1:nSites_spk,1,:)) .^ 2;
    vrVp = sum(mrVp);

    miSites_spk = single(P.miSites(1:nSites_spk, S0.spikeSites));
    mrX_spk = reshape(P.mrSiteXY(miSites_spk,1), size(miSites_spk));
    mrY_spk = reshape(P.mrSiteXY(miSites_spk,2), size(miSites_spk));

    mrPos_spk = zeros(numel(S0.spikeSites), 2, 'single');
    mrPos_spk(:,1) = sum(mrVp .* mrX_spk) ./ vrVp;
    mrPos_spk(:,2) = sum(mrVp .* mrY_spk) ./ vrVp;
end %func
