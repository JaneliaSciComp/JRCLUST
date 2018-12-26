function rescaleFigProj(hFigProj, projScale, hClust)
    %RESCALEFIGPROJ
    hCfg = hClust.hCfg;
    
%     if isnumeric(projScale)
%         S_fig.maxAmp = projScale;
%     else
%         S_fig.maxAmp = change_amp_(projScale, S_fig.maxAmp);
%     end

    rescaleUpdate(hFigProj, projScale, hCfg);

    if strcmp(hCfg.dispFeature, 'vpp')
        xLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        yLabel = 'Site # (%0.0f \\muV_{min})';
    elseif ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        xLabel = sprintf('Site # (PC %d)', hCfg.pcPair(1));
        yLabel = sprintf('Site # (PC %d)', hCfg.pcPair(2));
    else
        xLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        yLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end

    hFigProj.xlabel(sprintf(xLabel, projScale));
    hFigProj.ylabel(sprintf(yLabel, projScale));
end

%% LOCAL FUNCTIONS
function rescaleUpdate(hFigProj, projScale, hCfg)
    hFigProj.rmPlot('hSelect'); % clear select polygon
    bounds = projScale*[0 1];

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
