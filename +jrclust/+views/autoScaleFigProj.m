function autoScaleFigProj(hFigProj, hClust, selected)
    %AUTOSCALEFIGPROJ Automatically scale features in FigProj
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

    if all(cellfun(@isempty, projData))
        projScale = hFigProj.figData.boundScale;
    else
        projScale = max(cellfun(@(x) quantile(abs(x(:)), autoscalePct), projData));
    end

    jrclust.views.rescaleFigProj(hFigProj, projScale, hCfg);
end