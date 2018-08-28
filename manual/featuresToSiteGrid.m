%--------------------------------------------------------------------------
function [newX, newY, viPlot, tr_dim] = featuresToSiteGrid(xvals, yvals, bounds, maxPair, vpp)
    if nargin < 4
        maxPair = [];
    end
    if nargin < 5
        vpp = 1; % default view
    end

    % remap yvals and xvals into [0 1]
    xvals = linmap(xvals', bounds, [0 1], 1);
    yvals = linmap(yvals', bounds, [0 1], 1);

    [nSpikes, nSites] = size(yvals);
    if isempty(maxPair)
        maxPair = nSites;
    end

    % spike features translated into site-site boxes
    [transX, transY] = deal(nan([nSpikes, nSites, nSites], 'single'));

    for jSite = 1:nSites
        jSiteY = yvals(:, jSite);
        yMask = jSiteY > 0 & jSiteY < 1; % get points away from the boundaries

        for iSite = 1:nSites
            if abs(iSite - jSite) > maxPair
                continue;
            end

            if vpp && jSite > iSite % above diagonal: min vs. min
                iSiteX = yvals(:, iSite);
            else % diagonal and below: min vs. max
                iSiteX = xvals(:, iSite);
            end

            xMask = iSiteX > 0 & iSiteX < 1; % get points away from the boundaries

            xyMask = (xMask & yMask);
            transX(xyMask, jSite, iSite) = iSiteX(xyMask) + iSite - 1;
            transY(xyMask, jSite, iSite) = jSiteY(xyMask) + jSite - 1;
        end
    end

    % plot projection
    viPlot = find(~isnan(transX) & ~isnan(transY));
    newX = transX(viPlot);
    newX = newX(:);

    newY = transY(viPlot);
    newY = newY(:);
    tr_dim = size(transX);
end %func
