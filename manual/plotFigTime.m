%--------------------------------------------------------------------------
function plotFigTime(S0)
    % plot FigTime window. Uses subsampled data

    if nargin < 1
        S0 = get(0, 'UserData');
    end

    S_clu = S0.S_clu;
    P = S0.P;

    [hFig, S_fig] = getCachedFig('FigTime');

    %----------------
    % collect info
    iSite = S_clu.clusterSites(S0.iCluCopy);
    [backgroundFeatures, backgroundTimes] = getFet_site_(iSite, [], S0); % plot background
    [primaryFeatures, primaryTimes, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy, S0); % plot iCluCopy

    vcTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';
    if ~isempty(S0.iCluPaste)
        [secondaryFeatures, secondaryTimes] = getFet_site_(iSite, S0.iCluPaste, S0);
        vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', S0.iCluCopy, S0.iCluPaste, vcTitle);
    else
        secondaryFeatures = [];
        secondaryTimes = [];
        vcTitle = sprintf('Clu%d (black); %s', S0.iCluCopy, vcTitle);
    end
    time_lim = double([0, abs(S0.spikeTimes(end))] / P.sampleRateHz);

    %------------
    % draw
    if isempty(S_fig)
        S_fig.maxAmp = P.maxAmp;
        S_fig.hAx = axes_new_(hFig);
        set(S_fig.hAx, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % first time
        S_fig.hPlot0 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(1,:), 'MarkerSize', 5, 'LineStyle', 'none');
        S_fig.hPlot1 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
        S_fig.hPlot2 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(3,:), 'MarkerSize', 5, 'LineStyle', 'none'); %place holder
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
        S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
        if ~isempty(P.time_tick_show) %tick mark
            set(S_fig.hAx, 'XTick', time_lim(1):P.time_tick_show:time_lim(end));
        end
    end
    vpp_lim = [0, abs(S_fig.maxAmp)];
    % iFet = S_fig.iFet;
    % iFet = 1;
    if ~isfield(S_fig, 'iSite'), S_fig.iSite = []; end
    update_plot_(S_fig.hPlot0, backgroundTimes, backgroundFeatures);
    update_plot_(S_fig.hPlot1, primaryTimes, primaryFeatures);
    update_plot_(S_fig.hPlot2, secondaryTimes, secondaryFeatures);
    imrect_set_(S_fig.hRect, time_lim, vpp_lim);
    mouse_figure(hFig, S_fig.hAx); % allow zoom using wheel
    % button click function to select individual spikes, all spikes plotted

    if isfield(S_fig, 'vhAx_track')
        toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
        toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot1, S_fig.hPlot2, S_fig.hPlot0}, 1);
    end

    if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
    toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);

    axis_(S_fig.hAx, [time_lim, vpp_lim]);
    title_(S_fig.hAx, vcTitle);
    ylabel(S_fig.hAx, vcYlabel);

    S_fig = mergeStructs(S_fig, makeStruct_(iSite, time_lim, P, vpp_lim, viSpk1));
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
end %func
