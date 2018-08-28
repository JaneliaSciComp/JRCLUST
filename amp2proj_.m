%--------------------------------------------------------------------------
function [vrX, vrY, viPlot, tr_dim] = amp2proj_(yvals, xvals, bounds, maxPair)
    if nargin < 4
        maxPair = [];
    end

    % remap yvals and xvals into [0 1]
    xvals = linmap(xvals', bounds, [0 1], 1);
    yvals = linmap(yvals', bounds, [0 1], 1);

    [nSpikes, nSites] = size(yvals);
    if isempty(maxPair)
        maxPair = nSites;
    end

    [trX, trY] = deal(nan([nSpikes, nSites, nSites], 'single'));

    for jSite = 1:nSites
        vrY1 = yvals(:, jSite);
        yMask = vrY1 > 0 & vrY1 < 1; % get points away from the boundaries

        for iSite = 1:nSites
            if abs(iSite - jSite) > maxPair
                continue;
            end

            if jSite > iSite
                vrX1 = yvals(:, iSite);
            else
                vrX1 = xvals(:, iSite);
            end

            xMask = vrX1 > 0 & vrX1 < 1; % get points away from the boundaries

            xyMask = find(xMask & yMask);
            trX(xyMask, jSite, iSite) = vrX1(xyMask) + iSite - 1;
            trY(xyMask, jSite, iSite) = vrY1(xyMask) + jSite - 1;
        end
    end

    % plot projection
    viPlot = find(~isnan(trX) & ~isnan(trY));
    vrX = trX(viPlot);
    vrX = vrX(:);

    vrY = trY(viPlot);
    vrY = vrY(:);
    tr_dim = size(trX);
end %func
