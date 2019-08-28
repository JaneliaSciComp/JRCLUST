function [XData, YData, subset] = ampToProj(YData, XData, bounds, maxPair, hCfg)
    %AMPTOPROJ Reshape feature data from nSpikes x nSites to display in an
    %nSites x nSites grid
    %   input: YData, nSpikes x nSites, y-values for feature projection
    %   input: XData, nSpikes x nSites, x-values for feature projection
    [nSites, nSpikes] = size(YData);

    % subset features
    if nargout > 2
        subset = jrclust.utils.subsample(1:nSpikes, hCfg.nSpikesFigProj);
        YData = YData(:, subset);
        XData = XData(:, subset);
        nSpikes = numel(subset);
    end

    XData = linmap(XData', bounds, [0, 1]);
    YData = linmap(YData', bounds, [0, 1]);

    if isempty(maxPair)
        maxPair = nSites;
    end

    [jSites, iSites] = meshgrid(1:nSites, 1:nSites);
    maxPairMask = abs(jSites - iSites) <= maxPair;
    jSites = jSites(maxPairMask);
    iSites = iSites(maxPairMask);

    if strcmp(hCfg.dispFeature, 'vpp')
        aboveDiag = jSites > iSites;

        iSitesX = zeros(nSpikes, numel(iSites), 'like', XData);
        iSitesX(:, aboveDiag) = YData(:, iSites(aboveDiag));
        iSitesX(:, ~aboveDiag) = XData(:, iSites(~aboveDiag));
    else
        iSitesX = XData(:, iSites);
    end

    jSitesY = YData(:, jSites);

    ijMask = (jSitesY > 0 & jSitesY < 1) & (iSitesX > 0 & iSitesX < 1);

    iSitesX = iSitesX + (iSites - 1)';
    jSitesY = jSitesY + (jSites - 1)';

    XData = iSitesX(ijMask);
    XData = XData(:);

    YData = jSitesY(ijMask);
    YData = YData(:);

    if nargout > 2
        subset = repmat(subset', 1, size(ijMask, 2));
        subset = subset(ijMask);
    end
end

%% LOCAL FUNCTIONS
function vals = linmap(vals, oldLim, newLim)
    %LINMAP Rescale vals occurring within oldLim to newLim, saturating at
    %the boundaries
    if numel(oldLim) == 1
        oldLim = abs(oldLim)*[-1, 1];
    end

    % saturate at the boundaries
    vals(vals > oldLim(2)) = oldLim(2);
    vals(vals < oldLim(1)) = oldLim(1);

    if all(oldLim == 0) % nothing to rescale, all is 0 now
        return;
    end

    if oldLim(1) == oldLim(2) % ignore newLim and just rescale 
        vals = vals / oldLim(1);
    else
        vals = interp1(oldLim, newLim, vals, 'linear', 'extrap');
    end
end
