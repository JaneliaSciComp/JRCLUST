function rescaleFigProj(hFigProj, projScale, hCfg)
    %RESCALEFIGPROJ
    rescaleUpdate(hFigProj, projScale, hCfg);
    hFigProj.figData.boundScale = projScale;

    if strcmp(hCfg.dispFeature, 'vpp')
        XLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        YLabel = 'Site # (%0.0f \\muV_{min})';
    elseif ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        XLabel = sprintf('Site # (PC %d)', hCfg.pcPair(1));
        YLabel = sprintf('Site # (PC %d)', hCfg.pcPair(2));
    else
        XLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        YLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end

    hFigProj.axApply('default', @xlabel, sprintf(XLabel, projScale));
    hFigProj.axApply('default', @ylabel, sprintf(YLabel, projScale));
end

%% LOCAL FUNCTIONS
function rescaleUpdate(hFigProj, projScale, hCfg)
    if strcmp(hCfg.dispFeature, 'vpp')
        bounds = projScale*[0 1];
    else
        bounds = projScale*[-1 1];
    end

    dispFeatures = hFigProj.figData.dispFeatures;
    % rescale background features
    bgXData = dispFeatures.bgXData;
    bgYData = dispFeatures.bgYData;
    [bgXData, bgYData] = ampToProj(bgYData, bgXData, bounds, hCfg.nSiteDir, hCfg);
    hFigProj.updatePlot('background', bgXData, bgYData);

    % rescale foreground features
    fgYData = dispFeatures.fgYData;
    fgXData = dispFeatures.fgXData;
    [fgXData, fgYData] = ampToProj(fgYData, fgXData, bounds, hCfg.nSiteDir, hCfg);
    hFigProj.updatePlot('foreground', fgXData, fgYData);

    % rescale secondary foreground features
    fg2YData = dispFeatures.fg2YData;
    fg2XData = dispFeatures.fg2XData;
    if ~all(isnan([fg2XData(:); fg2YData(:)]))
        [fg2XData, fg2YData] = ampToProj(fg2YData, fg2XData, bounds, hCfg.nSiteDir, hCfg);
        hFigProj.updatePlot('foreground2', fg2XData, fg2YData);
    end
end
