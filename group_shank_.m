%--------------------------------------------------------------------------
function [spikeSites] = group_shank_(spikeSites, P)
    nSites = numel(P.chanMap);
    site2site = zeros([nSites, 1], 'like', spikeSites);
    [a,b,c] = unique(P.viShank_site);
    site2site(P.chanMap) = b(c);
    spikeSites = site2site(spikeSites);
end %func
