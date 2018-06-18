%--------------------------------------------------------------------------
function [viSite_spk] = group_shank_(viSite_spk, P)
    nSites = numel(P.viSite2Chan);
    site2site = zeros([nSites, 1], 'like', viSite_spk);
    [a,b,c] = unique(P.viShank_site);
    site2site(P.viSite2Chan) = b(c);
    viSite_spk = site2site(viSite_spk);
end %func
