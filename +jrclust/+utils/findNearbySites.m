function nearbySites = findNearbySites(siteLoc, iSite, evtDetectRad)
    %FINDNEARBYSITES Get sites which are within evtDetectRad of iSite
    dists = pdist2(siteLoc(iSite, :), siteLoc);
    nearbySites = find(dists <= evtDetectRad);
end
