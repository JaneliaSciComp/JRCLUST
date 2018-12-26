function autoScaleProjTime(hClust, hFigProj, hFigTime, selected)
    %AUTOSCALEPROJTIME Automatically scale features in FigProj and FigTime
    hCfg = hClust.hCfg;
    autoscalePct = hCfg.getOr('autoscalePct', 99.5)/100;

    % get display features for figProj, compute the autoscalePct'th
    % quantile, set boundaries of each site box to that scale
    dispFeaturesProj = hFigProj.figData.dispFeatures;
    fgYData = dispFeaturesProj.fgYData;
    fgXData = dispFeaturesProj.fgXData;
    fg2YData = dispFeaturesProj.fg2YData;
    fg2XData = dispFeaturesProj.fg2XData;

    if numel(selected) == 1
        projData = {fgYData, fgXData};
    else
        projData = {fgYData, fgXData, fg2YData, fg2XData};
    end
    projScale = max(cellfun(@(x) quantile(x(:), autoscalePct), projData));
    rescaleFigProj(hFigProj, projScale, hClust);

    % get display features for figTime on iSite, compute the
    % autoscalePct'th quantile, set ylimits to that scale
    iSite = hClust.clusterSites(selected(1));
    dispFeaturesTime = getFigTimeFeatures(hClust, iSite, selected(1));
    if numel(selected) == 1
        timeData = {dispFeaturesTime};
    else
        dispFeaturesTime2 = getFigTimeFeatures(hClust, iSite, selected(2));
        timeData = {dispFeaturesTime, dispFeaturesTime2};
    end
    timeScale = max(cellfun(@(x) quantile(x(:), autoscalePct), timeData));
    hFigTime.axSet('YLim', [0, 1]*timeScale);
    imrect_set_(hFigTime, 'hRect', [], [0, timeScale]);
end
