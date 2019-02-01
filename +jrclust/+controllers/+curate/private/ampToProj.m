function [XData, YData, validXY, boxShape] = ampToProj(YData, XData, bounds, maxPair, hCfg)
    %AMPTOPROJ Reshape feature data from nSpikes x nSites to display in an
    %nSites x nSites grid
    %   input: YData, nSpikes x nSites, y-values for feature projection
    %   input: XData, nSpikes x nSites, x-values for feature projection
    XData = linmap(XData', bounds, [0, 1]);
    YData = linmap(YData', bounds, [0, 1]);

    [nSpikes, nSites] = size(YData);
    if isempty(maxPair)
        maxPair = nSites;
    end

    % spike features translated into site-site boxes
    [boxedX, boxedY] = deal(nan([nSpikes, nSites, nSites], 'single'));

    for jSite = 1:nSites
        jSiteY = YData(:, jSite);
        yMask = jSiteY > 0  & jSiteY < 1; % get points away from the boundaries

        for iSite = 1:nSites
            if abs(iSite - jSite) > maxPair
                continue;
            end
            
            % vpp only:
            % min on site j vs. min on site i above the diagonal
            if strcmp(hCfg.dispFeature, 'vpp') && jSite > iSite
                isiteX = YData(:, iSite);
            else % diagonal and below: min vs. max
                isiteX = XData(:, iSite);
            end

            % points away from both x and y boundaries
            xyMask = (isiteX > 0 & isiteX < 1) & yMask;

            boxedX(xyMask, jSite, iSite) = isiteX(xyMask) + iSite - 1;
            boxedY(xyMask, jSite, iSite) = jSiteY(xyMask) + jSite - 1;
        end
    end

    validXY = find(~isnan(boxedX) & ~isnan(boxedY));
    XData = boxedX(validXY);
    XData = XData(:);

    YData = boxedY(validXY);
    YData = YData(:);

    boxShape = [nSpikes, nSites, nSites];
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
