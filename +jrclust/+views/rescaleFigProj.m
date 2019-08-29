function rescaleFigProj(hFigProj, projScale, hCfg)
    %RESCALEFIGPROJ
    rescaleUpdate(hFigProj, projScale, hCfg);
    hFigProj.figData.boundScale = projScale;

    if strcmp(hCfg.dispFeature, 'vpp')
        XLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        YLabel = 'Site # (%0.0f \\muV_{min})';
    elseif ~ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        XLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        YLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);        
    end

    if ~ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        hFigProj.axApply('default', @xlabel, sprintf(XLabel, projScale));
        hFigProj.axApply('default', @ylabel, sprintf(YLabel, projScale));
    end
end

%% LOCAL FUNCTIONS
function rescaleUpdate(hFigProj, projScale, hCfg)
    % get ratio between plotted and stored values
    s = hFigProj.figData.initialScale/projScale;

    % rescale background features
    bgXData = hFigProj.figData.background.XData(:);
    xFloor = floor(bgXData);
    if strcmp(hCfg.dispFeature, 'vpp')
        bgXData = (bgXData - xFloor)*s;
    else
        % remap to [-1, 1] first, then scale, then map back to [0, 1]
        bgXData = jrclust.utils.linmap(bgXData - xFloor, [0, 1], [-1, 1])*s;
        bgXData = jrclust.utils.linmap(bgXData, [-1, 1], [0, 1]);
    end
    bgXData((bgXData <= 0 | bgXData >= 1)) = nan;
    bgXData = bgXData + xFloor;

    bgYData = hFigProj.figData.background.YData(:);
    yFloor = floor(bgYData);
    if strcmp(hCfg.dispFeature, 'vpp')
        bgYData = (bgYData - yFloor)*s;
    else
        % remap to [-1, 1] first, then scale, then map back to [0, 1]
        bgYData = jrclust.utils.linmap(bgYData - yFloor, [0, 1], [-1, 1])*s;
        bgYData = jrclust.utils.linmap(bgYData, [-1, 1], [0, 1]);
    end
    bgYData((bgYData <= 0 | bgYData >= 1)) = nan;
    bgYData = bgYData + yFloor;
    hFigProj.updatePlot('background', bgXData, bgYData);

    % rescale foreground features
    subset = hFigProj.figData.foreground.subset; % just scale the plotted subset
	fgXData = hFigProj.figData.foreground.XData(subset, :);
    fgXData = fgXData(:);
	xFloor = floor(fgXData);
    if strcmp(hCfg.dispFeature, 'vpp')
        fgXData = (fgXData - xFloor)*s;
    else
        % remap to [-1, 1] first, then scale, then map back to [0, 1]
        fgXData = jrclust.utils.linmap(fgXData - xFloor, [0, 1], [-1, 1])*s;
        fgXData = jrclust.utils.linmap(fgXData, [-1, 1], [0, 1]);
    end
    fgXData((fgXData <= 0 | fgXData >= 1)) = nan;
    fgXData = fgXData + xFloor;

    fgYData = hFigProj.figData.foreground.YData(subset, :);
    fgYData = fgYData(:);
    yFloor = floor(fgYData);
    if strcmp(hCfg.dispFeature, 'vpp')
        fgYData = (fgYData - yFloor)*s;
    else
        % remap to [-1, 1] first, then scale, then map back to [0, 1]
        fgYData = jrclust.utils.linmap(fgYData - yFloor, [0, 1], [-1, 1])*s;
        fgYData = jrclust.utils.linmap(fgYData, [-1, 1], [0, 1]);
    end
    fgYData((fgYData <= 0 | fgYData >= 1)) = nan;
    fgYData = fgYData + yFloor;
    hFigProj.updatePlot('foreground', fgXData, fgYData);

    % rescale secondary foreground features
    fg2XData = hFigProj.figData.foreground2.XData(:);
	fg2YData = hFigProj.figData.foreground2.YData(:);

    if ~all(isnan([fg2XData(:); fg2YData(:)]))
        xFloor = floor(fg2XData);
        if strcmp(hCfg.dispFeature, 'vpp')
            fg2XData = (fg2XData - xFloor)*s;
        else
            % remap to [-1, 1] first, then scale, then map back to [0, 1]
            fg2XData = jrclust.utils.linmap(fg2XData - xFloor, [0, 1], [-1, 1])*s;
            fg2XData = jrclust.utils.linmap(fg2XData, [-1, 1], [0, 1]);
        end
        fg2XData((fg2XData <= 0 | fg2XData >= 1)) = nan;
        fg2XData = fg2XData + xFloor;

        yFloor = floor(fg2YData);
        if strcmp(hCfg.dispFeature, 'vpp')
            fg2YData = (fg2YData - yFloor)*s;
        else
            % remap to [-1, 1] first, then scale, then map back to [0, 1]
            fg2YData = jrclust.utils.linmap(fg2YData - yFloor, [0, 1], [-1, 1])*s;
            fg2YData = jrclust.utils.linmap(fg2YData, [-1, 1], [0, 1]);
        end
        fg2YData((fg2YData <= 0 | fg2YData >= 1)) = nan;
        fg2YData = fg2YData + yFloor;
        hFigProj.updatePlot('foreground2', fg2XData, fg2YData);
    end
end
