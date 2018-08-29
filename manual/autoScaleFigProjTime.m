%--------------------------------------------------------------------------
function autoScaleFigProjTime(S0, fPlot)
    % auto-scale and refresh

    if nargin < 1
        S0 = get(0, 'UserData');
    end
    if nargin < 2
        fPlot = 0;
    end

    autoscale_pct = getOr(S0.P, 'autoscale_pct', 99.5);

    [hFigProj, figProjData] = getCachedFig('FigProj');
    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = getFigProjFeatures(S0, figProjData.sitesOfInterest);

    if isempty(mrMin2) || isempty(mrMax2)
        featureData = {mrMin1, mrMax1};
    else
        featureData = {mrMin1, mrMax1, mrMin2, mrMax2};
    end

    % set maxAmp to max of autoscale_pct percentile, taken over all relevant coordinates
    figProjData.maxAmp = max(cellfun(@(x) quantile(abs(x(:)), autoscale_pct/100), featureData));
    set(hFigProj, 'UserData', figProjData);

    % Update time
    [hFigTime, figTimeData] = getCachedFig('FigTime');
    primaryClusterSite = S0.S_clu.clusterSites(S0.primarySelectedCluster);

    [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(primaryClusterSite, S0.primarySelectedCluster, S0); % plot primarySelectedCluster
    if isempty(S0.secondarySelectedCluster)
        cvrFet = {vrFet1};
    else
        [vrFet2, vrTime2, vcYlabel, viSpk2] = getFet_site_(primaryClusterSite, S0.secondarySelectedCluster, S0); % plot primarySelectedCluster
        cvrFet = {vrFet1, vrFet2};
    end

    % set maxAmp to max of autoscale_pct percentile, taken over all foreground amplitudes on site
    figTimeData.maxAmp = max(cellfun(@(x) quantile(x(:), autoscale_pct/100), cvrFet));
    set(hFigTime, 'UserData', figTimeData);

    % plot
    if fPlot
        keyPressFcn_cell_(getCachedFig('FigWav'), {'j', 't'}, S0);
    else
        rescaleFigProj(figProjData.maxAmp, hFigProj, figProjData, S0);
        rescale_FigTime_(figTimeData.maxAmp, S0, S0.P);
    end
end % function
