function nearbySites = findNearSite_(siteLoc, iSite, evtDetectRad)
    dists = pdist2_(siteLoc(iSite, :), siteLoc);
    nearbySites = find(dists <= evtDetectRad);
end
