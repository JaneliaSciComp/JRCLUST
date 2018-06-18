%--------------------------------------------------------------------------
function [mrPos_spk, viSpk_re] = position_spk_(viSite_spk, tnWav_spk, P)

    mrPos_site1 = P.mrSiteXY(P.miSites(1, viSite_spk), :)'; %first pos

    % determine centroid location and second largest amplitude
    [mrPos_spk, mrA_spk] = tnWav_centroid_(tnWav_spk, viSite_spk, P);
    [~, viiSite2_spk] = max(mrA_spk((2:end), :)); %find second max
    miSites2 = P.miSites(2:end, :);
    viiSite2_spk = sub2ind(size(miSites2), viiSite2_spk(:), viSite_spk(:));
    viSite_spk2 = miSites2(viiSite2_spk);

    % Find where second largest site is closer to the spike centroid
    mrPos_site2 = P.mrSiteXY(viSite_spk2, :)';
    dist__ = @(mr1,mr2)sum((mr1-mr2).^2);
    viSpk_re = find(dist__(mrPos_spk,mrPos_site2) < dist__(mrPos_spk, mrPos_site1));
end %func
