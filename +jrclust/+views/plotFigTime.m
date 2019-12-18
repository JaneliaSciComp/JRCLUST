function hFigTime = plotFigTime(hFigTime, hClust, hCfg, selected, maxAmp, iSite)
    %DOPLOTFIGTIME Plot features vs. time
    timeLimits = double([0, abs(hClust.spikeTimes(end))/hCfg.sampleRate]);

    % construct plot for the first time
    if ~hFigTime.hasAxes('default')
        hFigTime.addAxes('default');
        hFigTime.axApply('default', @set, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % first time
        hFigTime.addPlot('background', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(1, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addPlot('foreground', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(2, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addPlot('foreground2', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(3, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.axApply('default', @xlabel, 'Time (s)');
        hFigTime.axApply('default', @grid, 'on');

        % rectangle plot
        rectPos = [timeLimits(1), maxAmp, diff(timeLimits), maxAmp];
        hFigTime.addPlot('hRect', @imrect, rectPos);
        hFigTime.plotApply('hRect', @setColor, 'r');
        hFigTime.plotApply('hRect', @setPositionConstraintFcn, makeConstrainToRectFcn('imrect', timeLimits, [-4000 4000]));

        hFigTime.setHideOnDrag('background'); % hide background spikes when dragging
    end

    [bgFeatures, bgTimes] = getFigTimeFeatures(hClust, iSite); % plot background
    [fgFeatures, fgTimes, YLabel] = getFigTimeFeatures(hClust, iSite, selected(1)); % plot primary selected cluster

    if numel(selected) == 2
        [fgFeatures2, fgTimes2] = getFigTimeFeatures(hClust, iSite, selected(2));
        figTitle = sprintf('Unit %d (black), Unit %d (red); (press [H] for help)', selected(1), selected(2));
    else
        fgFeatures2 = [];
        fgTimes2 = [];
        figTitle = sprintf('Unit %d (black); (press [H] for help)', selected(1));
    end

    vppLim = [0, abs(maxAmp)];

    hFigTime.updatePlot('background', bgTimes, bgFeatures);
    hFigTime.updatePlot('foreground', fgTimes, fgFeatures);
    hFigTime.updatePlot('foreground2', fgTimes2, fgFeatures2);
    imrectSetPosition(hFigTime, 'hRect', timeLimits, vppLim);

%     if isfield(S_fig, 'vhAx_track')
%         toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
%         toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot1, S_fig.hPlot2, S_fig.hPlot0}, 1);
%     end
% 
    if ~isfield(hFigTime.figData, 'doPlotBG')
        hFigTime.figData.doPlotBG = 1;
    end

    hFigTime.axApply('default', @axis, [timeLimits, vppLim]);
    hFigTime.axApply('default', @title, figTitle);
    hFigTime.axApply('default', @ylabel, YLabel);
end