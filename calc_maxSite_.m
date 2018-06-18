%--------------------------------------------------------------------------
% 10/25/17 JJJ: Created and tested
function P = calc_maxSite_(P)
    % Auto determine maxSite from the radius and site ifno

    P.maxSite = get_(P, 'maxSite');
    P.nSites_ref = get_(P, 'nSites_ref');
    if ~isempty(P.maxSite) && ~isempty(P.nSites_ref), return; end
    mrDist_site = pdist2(P.mrSiteXY, P.mrSiteXY);
    maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
    nSites_fet = max(sum(mrDist_site <= maxDist_site_um)); % 11/7/17 JJJ: med to max
    if isempty(P.nSites_ref)
        maxDist_site_spk_um = get_set_(P, 'maxDist_site_spk_um', maxDist_site_um+25);
        nSites_spk = max(sum(mrDist_site <= maxDist_site_spk_um)); % 11/7/17 JJJ: med to max
        maxSite = (nSites_spk-1)/2;
        P.nSites_ref = nSites_spk - nSites_fet;
    else
        nSites_spk = nSites_fet + P.nSites_ref;
        maxSite = (nSites_spk-1)/2;
    end

    if isempty(P.maxSite), P.maxSite = maxSite; end
    fprintf('Auto-set: maxSite=%0.1f, nSites_ref=%0.1f\n', P.maxSite, P.nSites_ref);
end %func
