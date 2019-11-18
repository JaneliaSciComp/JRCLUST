function autoScaleFigTime(hFigTime, hClust, selected, iSite, channel_idx)
    %AUTOSCALEFIGTIME Automatically scale features in FigTime
    hCfg = hClust.hCfg;
    autoscalePct = hCfg.getOr('autoscalePct', 99.5)/100;

    % get display features for figTime on iSite, compute the
    % autoscalePct'th quantile, set ylimits to that scale
    dispFeaturesTime = getFigTimeFeatures(hClust, iSite, selected(1),channel_idx);
    if numel(dispFeaturesTime)<2 || ~any(dispFeaturesTime)
        dispFeaturesTime = getFigTimeFeatures(hClust, iSite,[],channel_idx);
    end
    if numel(selected) == 1
        timeData = {dispFeaturesTime};
    else
        dispFeaturesTime2 = getFigTimeFeatures(hClust, iSite, selected(2),channel_idx);
        timeData = {dispFeaturesTime, dispFeaturesTime2};
    end

    if all(cellfun(@isempty, timeData))
        return;
    end

    timeScale = max(cellfun(@(x) quantile(x(:), autoscalePct), timeData));
    if isnan(timeScale)
        timeScale = 1;
    end
    hFigTime.axApply('default', @set, 'YLim', [0, 1]*timeScale);
    hFigTime.axApply('histogram',@set,'YLim',[0,1]*timeScale);
    imrectSetPosition(hFigTime, 'hRect', [], [0, timeScale]);
end