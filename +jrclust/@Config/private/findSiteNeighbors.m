function [sites, dists] = findSiteNeighbors(siteLoc, nNeighbors, ignoreSites, shankMap)
    %FINDSITENEIGHBORS Get the nearest neighbors for each site
    if nargin < 3
        ignoreSites = [];
    end
    if nargin < 4 || numel(unique(shankMap)) <= 1
        shankMap = [];
    end

    % set ignored sites infinitely far away
    if ~isempty(ignoreSites)
        siteLoc(ignoreSites, :) = inf;
    end

    nSites = size(siteLoc, 1);
    nNeighbors = min(nNeighbors, nSites);

    [sites, dists] = deal(zeros(nNeighbors, nSites));

    % matrix of pairwise distances of site locations
    pdists = pdist2(siteLoc, siteLoc);
    for iSite = 1:nSites
        siteDists = pdists(iSite, :);

        % set sites on a different shank infinitely far away
        if ~isempty(shankMap)
            siteDists(shankMap ~= shankMap(iSite)) = inf;
        end

        [siteDists, iSiteDists] = sort(siteDists, 'ascend');

        % take the nNeighbors nearest sites, in order
        sites(:, iSite) = iSiteDists(1:nNeighbors);
        dists(:, iSite) = siteDists(1:nNeighbors);
    end
end
