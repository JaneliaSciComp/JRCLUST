%--------------------------------------------------------------------------
function viSiteNear = findNearSite_(mrSiteXY, iSite, maxDist_site_um)
    vrDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
    viSiteNear = find(vrDist <= maxDist_site_um);
end %func
