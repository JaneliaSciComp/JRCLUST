function [XData, YData] = ampToProj(YData, XData, bounds, maxPair, hCfg)
    %AMPTOPROJ Reshape feature data from nSpikes x nSites to display in an
    %nSites x nSites grid
    %   input: YData, nSites x nSpikes, y-values for feature projection
    %   input: XData, nSites x nSpikes, x-values for feature projection
    [nSites, nSpikes] = size(YData);

    XData = jrclust.utils.linmap(XData', bounds, [0, 1]);
    YData = jrclust.utils.linmap(YData', bounds, [0, 1]);

    if isempty(maxPair)
        maxPair = nSites;
    end

    [jSites, iSites] = meshgrid(1:nSites, 1:nSites);
    maxPairMask = abs(jSites - iSites) <= maxPair;
    jSites = jSites(maxPairMask);
    iSites = iSites(maxPairMask);

    % when using VPP as a feature, above diagonal is vmin vs. vmin, whereas
    % on the diagonal and below is vmin vs. vmax
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

    iSitesX(~ijMask) = nan;
    XData = iSitesX;

    jSitesY(~ijMask) = nan;
    YData = jSitesY;
end
