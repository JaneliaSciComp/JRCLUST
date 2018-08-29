%--------------------------------------------------------------------------
function plotFigTime(S0)
    % plot FigTime window. Uses subsampled data

    if nargin<1, S0 = get(0, 'UserData'); end
    S_clu = S0.S_clu; P = S0.P;
    [hFig, S_fig] = getCachedFig('FigTime');

    %----------------
    % collect info
    iSite = S_clu.clusterSites(S0.primarySelectedCluster);
    [vrFet0, vrTime0] = getFet_site_(iSite, [], S0); % plot background
    [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.primarySelectedCluster, S0); % plot primarySelectedCluster

    vcTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';
    if ~isempty(S0.secondarySelectedCluster)
        [vrFet2, vrTime2] = getFet_site_(iSite, S0.secondarySelectedCluster, S0);
        vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', S0.primarySelectedCluster, S0.secondarySelectedCluster, vcTitle);
    else
        vrFet2 = [];
        vrTime2 = [];
        vcTitle = sprintf('Clu%d (black); %s', S0.primarySelectedCluster, vcTitle);
    end
    time_lim = double([0, abs(S0.spikeTimes(end))] / P.sampleRateHz);

    %------------
    % draw
    if isempty(S_fig)
        S_fig.maxAmp = P.maxAmp;
        S_fig.hAx = newAxes(hFig);
        set(S_fig.hAx, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % first time
        S_fig.hPlotBG = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(1,:), 'MarkerSize', 5, 'LineStyle', 'none');
        S_fig.hPlotFG = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
        S_fig.hPlotFG2 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(3,:), 'MarkerSize', 5, 'LineStyle', 'none'); %place holder
        xlabel('Time (s)');
        grid on;

        % rectangle plot
        vrPos_rect = [time_lim(1), S_fig.maxAmp, diff(time_lim), S_fig.maxAmp];
        S_fig.hRect = imrect_(S_fig.hAx, vrPos_rect); %default position?
        if ~isempty(S_fig.hRect)
            setColor(S_fig.hRect, 'r');
            setPositionConstraintFcn(S_fig.hRect, ...
            makeConstrainToRectFcn('imrect',time_lim, [-4000 4000]));
        end
        set(hFig, 'KeyPressFcn', @keyPressFigTime);
        S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlotBG, S_fig);
        if ~isempty(P.time_tick_show) %tick mark
            set(S_fig.hAx, 'XTick', time_lim(1):P.time_tick_show:time_lim(end));
        end
    end
    vpp_lim = [0, abs(S_fig.maxAmp)];
    % iFet = S_fig.iFet;
    % iFet = 1;
    if ~isfield(S_fig, 'iSite'), S_fig.iSite = []; end
    updatePlot(S_fig.hPlotBG, vrTime0, vrFet0);
    updatePlot(S_fig.hPlotFG, vrTime1, vrFet1);
    updatePlot(S_fig.hPlotFG2, vrTime2, vrFet2);
    imrect_set_(S_fig.hRect, time_lim, vpp_lim);
    mouse_figure(hFig, S_fig.hAx); % allow zoom using wheel
    % button click function to select individual spikes, all spikes plotted

    if isfield(S_fig, 'vhAx_track')
        toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
        toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlotFG, S_fig.hPlotFG2, S_fig.hPlotBG}, 1);
    end

    if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
    toggleVisible_(S_fig.hPlotBG, S_fig.fPlot0);

    axis_(S_fig.hAx, [time_lim, vpp_lim]);
    title_(S_fig.hAx, vcTitle);
    ylabel(S_fig.hAx, vcYlabel);

    S_fig = mergeStructs(S_fig, makeStruct(iSite, time_lim, P, vpp_lim, viSpk1));
    S_fig.csHelp = {...
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

    set(hFig, 'UserData', S_fig);
end % function
