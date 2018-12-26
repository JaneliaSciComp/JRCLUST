function hFigTime = doPlotFigTime(hFigTime, hClust, hCfg, selected, maxAmp, iSite)
    %DOPLOTFIGTIME
    if nargin < 6
        iSite = hClust.clusterSites(selected(1));
    end

    timeLimits = double([0, abs(hClust.spikeTimes(end))/hCfg.sampleRate]);
    % construct plot for the first time
    if isempty(hFigTime.figData)
        hFigTime.axes();
        hFigTime.axSet('Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % first time
        hFigTime.addLine('background',  nan, nan, 'Marker', '.', 'Color', hCfg.mrColor_proj(1, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addLine('foreground',  nan, nan, 'Marker', '.', 'Color', hCfg.mrColor_proj(2, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addLine('foreground2', nan, nan, 'Marker', '.', 'Color', hCfg.mrColor_proj(3, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.xlabel('Time (s)');
        hFigTime.grid('on');

        % rectangle plot
        rectPos = [timeLimits(1), maxAmp, diff(timeLimits), maxAmp];
        hFigTime.addImrect('hRect', rectPos);
        hFigTime.imrectFun('hRect', @setColor, 'r');
        hFigTime.imrectFun('hRect', @setPositionConstraintFcn, makeConstrainToRectFcn('imrect', timeLimits, [-4000 4000]));

        hFigTime.setHideOnDrag('background'); % hide background spikes when dragging
        if ~isempty(hCfg.time_tick_show) % tick mark
            hFigTime.axSet('XTick', timeLimits(1):hCfg.time_tick_show:timeLimits(end));
        end
    end

    [bgFeatures, bgTimes] = getFigTimeFeatures(hClust, iSite); % plot background
    [fgFeatures, fgTimes, yLabel] = getFigTimeFeatures(hClust, iSite, selected(1)); % plot primary selected cluster

    figTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';
    if numel(selected) == 2
        [fgFeatures2, fgTimes2] = getFigTimeFeatures(hClust, iSite, selected(2));
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', selected(1), selected(2), figTitle);
    else
        fgFeatures2 = [];
        fgTimes2 = [];
        figTitle = sprintf('Clu%d (black); %s', selected(1), figTitle);
    end

    vppLim = [0, abs(maxAmp)];

    hFigTime.updatePlot('background', bgTimes, bgFeatures);
    hFigTime.updatePlot('foreground', fgTimes, fgFeatures);
    hFigTime.updatePlot('foreground2', fgTimes2, fgFeatures2);
    imrect_set_(hFigTime, 'hRect', timeLimits, vppLim);

%     if isfield(S_fig, 'vhAx_track')
%         toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
%         toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot1, S_fig.hPlot2, S_fig.hPlot0}, 1);
%     end
% 
    if ~isfield(hFigTime.figData, 'doPlotBG')
        hFigTime.figData.doPlotBG = true;
    end
%     toggleVisible_(S_fig.hPlot0, S_fig.doPlotBG);

    hFigTime.axis([timeLimits, vppLim]);
    hFigTime.title(figTitle);
    hFigTime.ylabel(yLabel);

%     S_fig = struct_merge_(S_fig, makeStruct_(iSite, timeLimits, hCfg, vpp_lim, clusterSpikes));
    hFigTime.figData.csHelp = {'Up/Down: change channel', ...
                               'Left/Right: Change sites', ...
                               'Shift + Left/Right: Show different features', ...
                               'r: reset scale', ...
                               'a: auto-scale', ...
                               'c: show pca across sites', ...
                               'e: export cluster info', ...
                               'f: export cluster feature', ...
                               'Zoom: mouse wheel', ...
                               'H-Zoom: press x and wheel. space to reset', ...
                               'V-Zoom: press y and wheel. space to reset', ...
                               'Drag while pressing wheel: pan'};
end
