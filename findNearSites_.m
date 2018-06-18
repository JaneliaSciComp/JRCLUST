%--------------------------------------------------------------------------
function [miSites, mrDist] = findNearSites_(mrSiteXY, maxSite, viSiteZero, viShank_site)
    % find nearest sites
    if nargin < 3, viSiteZero = []; end
    if nargin<4, viShank_site = []; end
    if numel(unique(viShank_site)) <= 1, viShank_site = []; end
    max_dist = max(pdist(mrSiteXY));
    if ~isempty(viSiteZero)
        mrSiteXY(viSiteZero,:) = max_dist*2; %bad sites will never be near
    end
    nNearSites = maxSite*2+1;
    nSites = size(mrSiteXY,1);
    nNearSites = min(nNearSites, nSites);
    [miSites, mrDist] = deal(zeros(nNearSites, nSites));
    viNearSites = 1:nNearSites;
    for iSite=1:nSites
        vrSiteDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
        if ~isempty(viShank_site)
            vrSiteDist(viShank_site(iSite) ~= viShank_site) = max_dist*4;
        end
        [vrSiteDist, viSrt] = sort(vrSiteDist, 'ascend');
        miSites(:,iSite) = viSrt(viNearSites);
        mrDist(:,iSite) = vrSiteDist(viNearSites);
    end
end %func
