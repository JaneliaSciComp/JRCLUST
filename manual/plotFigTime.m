%--------------------------------------------------------------------------
function plotFigTime(S0)
    % plot FigTime window. Uses subsampled data

    if nargin < 1
        S0 = get(0, 'UserData');
    end
    S_clu = S0.S_clu;
    P = S0.P;
    [hFig, figData] = getCachedFig('FigTime');

    %----------------
    % collect info
    primarySite = S_clu.clusterSites(S0.primarySelectedCluster);
    [backgroundFeatures, backgroundTimes] = getFigTimeFeatures(primarySite, [], S0); % plot background
    [foregroundFeatures, foregroundTimes, vcYlabel, viSpk1] = getFigTimeFeatures(primarySite, S0.primarySelectedCluster, S0); % plot primarySelectedCluster

    figTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';

    if ~isempty(S0.secondarySelectedCluster)
        [vrFet2, vrTime2] = getFigTimeFeatures(primarySite, S0.secondarySelectedCluster, S0);
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', S0.primarySelectedCluster, S0.secondarySelectedCluster, figTitle);
    else
        vrFet2 = [];
        vrTime2 = [];
        figTitle = sprintf('Clu%d (black); %s', S0.primarySelectedCluster, figTitle);
    end
    timeLimitSecs = double([0, abs(S0.spikeTimes(end))] / P.sampleRateHz);

    %------------
    % draw
    if isempty(figData)
        figData.maxAmp = P.maxAmp;
        figData.hAx = newAxes(hFig);
        set(figData.hAx, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % first time
        figData.hPlotBG = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(1,:), 'MarkerSize', 5, 'LineStyle', 'none');
        figData.hPlotFG = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
        figData.hPlotFG2 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(3,:), 'MarkerSize', 5, 'LineStyle', 'none'); %place holder
        xlabel('Time (s)');
        grid on;

        % rectangle plot
        vrPos_rect = [timeLimitSecs(1), figData.maxAmp, diff(timeLimitSecs), figData.maxAmp];
        figData.hRect = imrect_(figData.hAx, vrPos_rect); %default position?
        if ~isempty(figData.hRect)
            setColor(figData.hRect, 'r');
            setPositionConstraintFcn(figData.hRect, ...
            makeConstrainToRectFcn('imrect',timeLimitSecs, [-4000 4000]));
        end
        set(hFig, 'KeyPressFcn', @keyPressFigTime);
        figData.cvhHide_mouse = mouse_hide_(hFig, figData.hPlotBG, figData);
        if ~isempty(P.time_tick_show) %tick mark
            set(figData.hAx, 'XTick', timeLimitSecs(1):P.time_tick_show:timeLimitSecs(end));
        end
    end

    vpp_lim = [0, abs(figData.maxAmp)];

    if ~isfield(figData, 'primarySite')
        figData.primarySite = [];
    end

    updatePlot(figData.hPlotBG, backgroundTimes, backgroundFeatures);
    updatePlot(figData.hPlotFG, foregroundTimes, foregroundFeatures);
    updatePlot(figData.hPlotFG2, vrTime2, vrFet2);
    imrect_set_(figData.hRect, timeLimitSecs, vpp_lim);
    mouse_figure(hFig, figData.hAx); % allow zoom using wheel
    % button click function to select individual spikes, all spikes plotted

    if isfield(figData, 'vhAx_track')
        toggleVisible_({figData.vhAx_track, figData.hPlot0_track, figData.hPlot1_track, figData.hPlot2_track}, 0);
        toggleVisible_({figData.hAx, figData.hRect, figData.hPlotFG, figData.hPlotFG2, figData.hPlotBG}, 1);
    end

    if ~isfield(figData, 'fPlot0')
        figData.fPlot0 = 1;
    end
    toggleVisible_(figData.hPlotBG, figData.fPlot0);

    axis_(figData.hAx, [timeLimitSecs, vpp_lim]);
    title_(figData.hAx, figTitle);
    ylabel(figData.hAx, vcYlabel);

    figData = mergeStructs(figData, makeStruct(primarySite, timeLimitSecs, P, vpp_lim, viSpk1));
    figData.csHelp = {...
        'Up/Down: change channel', ...
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

    set(hFig, 'UserData', figData);
end % function
