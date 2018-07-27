%--------------------------------------------------------------------------
function nCols = nColumns_probe_(P)
    % Checkerboard four-column is considered as two column probe since
    % two sites per vertical step
    viShank_site = get_(P, 'viShank_site');
    if ~isempty(viShank_site)
        viSites = find(P.viShank_site == P.viShank_site(1));
        vrSiteY = P.mrSiteXY(viSites,2);
    else
        vrSiteY = P.mrSiteXY(:,2);
    end
    vrSiteY_unique = unique(vrSiteY);
    vnSites_group = hist(vrSiteY, vrSiteY_unique);
    nCols = median(vnSites_group);
end %func
