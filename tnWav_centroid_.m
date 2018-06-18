%--------------------------------------------------------------------------
function [mrPos_spk, mrVpp_spk] = tnWav_centroid_(tnWav_spk, viSite_spk, P)
    % 2 x nSpk matrix containing centroid and amplitude

    mrVpp_spk = trWav2fet_(tnWav_spk, P); % apply car
    % mrVpp_spk = tr2Vpp(tnWav_spk, P)';

    mrVpp_spk1 = mrVpp_spk .^ 2; % compute covariance with the center site
    mrVpp_spk1 = bsxfun(@rdivide, mrVpp_spk1, sum(mrVpp_spk1, 1));
    miSite_spk = P.miSites(:, viSite_spk);
    mrSiteXY = single(P.mrSiteXY);
    mrSiteX_spk = reshape(mrSiteXY(miSite_spk(:), 1), size(miSite_spk));
    mrSiteY_spk = reshape(mrSiteXY(miSite_spk(:), 2), size(miSite_spk));
    mrPos_spk = [sum(mrVpp_spk1 .* mrSiteX_spk, 1); sum(mrVpp_spk1 .* mrSiteY_spk, 1)];
end %func
