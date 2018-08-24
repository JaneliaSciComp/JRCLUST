%--------------------------------------------------------------------------
function [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, maxPair)
    if nargin < 4
        maxPair = [];
    end

    mrMax = linmap_(mrMax', [0, 1] * maxAmp, [0,1], 1);
    mrMin = linmap_(mrMin', [0, 1] * maxAmp, [0,1], 1);

    [nSpikes, nSites] = size(mrMin);
    if isempty(maxPair)
        maxPair = nSites;
    end

    [trX, trY] = deal(nan([nSpikes, nSites, nSites], 'single'));

    for jSite = 1:nSites
        vrY1 = mrMin(:, jSite);
        yMask = vrY1 > 0 & vrY1 < 1;

        for iSite = 1:nSites
            if abs(iSite - jSite) > maxPair
                continue;
            end

            if jSite > iSite
                vrX1 = mrMin(:, iSite);
            else
                vrX1 = mrMax(:, iSite);
            end

            viPlot1 = find(vrX1 > 0 & vrX1 < 1 & yMask);
            trX(viPlot1,jSite,iSite) = vrX1(viPlot1) + iSite - 1;
            trY(viPlot1,jSite,iSite) = vrY1(viPlot1) + jSite - 1;
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
