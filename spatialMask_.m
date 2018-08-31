%--------------------------------------------------------------------------
% 11/7/17 JJJ: Created
function [vrWeight_site1, vrDist_site1] = spatialMask_(P, iSite, nSites_spk, decayDist_um)
    if nargin<3, nSites_spk = size(P.miSites,1); end
    if nargin<4, decayDist_um = P.maxDist_site_um; end

    mrSiteXY1 = P.mrSiteXY(P.miSites(1:nSites_spk, iSite),:);
    vrDist_site1 = pdist2(mrSiteXY1(1,:), mrSiteXY1);
    vrWeight_site1 = 2.^(-vrDist_site1 / decayDist_um); %P.maxDist_site_merge_um);
    vrWeight_site1 = vrWeight_site1(:);
    vrDist_site1 = vrDist_site1(:);
end % function
