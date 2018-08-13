%--------------------------------------------------------------------------
function plotFigProj(S0)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    S_clu = S0.S_clu;
    P = S0.P;

    [hFig, S_fig] = getCachedFig('FigProj');

    primaryCluster = S0.iCluCopy;
    secondaryCluster = S0.iCluPaste;
    update_plot2_proj_(); % erase prev objects

    %---------------
    % Compute
    iSite1 = S_clu.clusterSites(primaryCluster);

    if ~isfield(P, 'sitesOfInterest')
        P.sitesOfInterest = sort(P.miSites(:, iSite1), 'ascend');
    end

    sitesOfInterest = P.sitesOfInterest;
    nSites = numel(P.sitesOfInterest);
    cell_plot = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};

    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
            xLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
            yLabel = 'Site # (%0.0f \\muV_{min})';
        otherwise
            xLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.displayFeature, P.displayFeature, P.displayFeature);
            yLabel = sprintf('Site # (%%0.0f %s)', P.displayFeature);
    end
    figTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

    %----------------
    % display
    if isempty(S_fig)
        S_fig.maxAmp = P.maxAmp;
        S_fig.hAx = axes_new_(hFig);
        set(S_fig.hAx, 'Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');

        S_fig.hPlot0 = line(nan, nan, 'Color', P.mrColor_proj(1,:), 'Parent', S_fig.hAx);
        S_fig.hPlot1 = line(nan, nan, 'Color', P.mrColor_proj(2,:), 'Parent', S_fig.hAx); %place holder
        S_fig.hPlot2 = line(nan, nan, 'Color', P.mrColor_proj(3,:), 'Parent', S_fig.hAx); %place holder
        set([S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2], cell_plot{:}); %common style

        S_fig.sitesOfInterest = []; % so that it can update
        S_fig.displayFeature = 'vpp';

        % plot boundary
        plotTable_([0, nSites], '-', 'Color', [.5 .5 .5]); %plot in one scoop
        plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5); %plot in one scoop

        mouse_figure(hFig);
        set(hFig, 'KeyPressFcn', @keyPressFcn_FigProj_);
        S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
        set_fig_(hFig, S_fig);
    end

    % get features for x0, y0, S_plot0 in one go
    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, P.sitesOfInterest);

    if ~isfield(S_fig, 'sitesOfInterest')
        S_fig.sitesOfInterest = [];
    end

    if ~equal_vr_(S_fig.sitesOfInterest, P.sitesOfInterest)
        plot_proj_(S_fig.hPlot0, mrMin0, mrMax0, P, S_fig.maxAmp);
    end

    plot_proj_(S_fig.hPlot1, mrMin1, mrMax1, P, S_fig.maxAmp);

    if ~isempty(secondaryCluster)
        plot_proj_(S_fig.hPlot2, mrMin2, mrMax2, P, S_fig.maxAmp);
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', primaryCluster, secondaryCluster, figTitle);
    else
        update_plot_(S_fig.hPlot2, nan, nan);
        figTitle = sprintf('Clu%d (black); %s', primaryCluster, figTitle);
    end

    % Annotate axes
    axis_(S_fig.hAx, [0 nSites 0 nSites]);
    set(S_fig.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', P.sitesOfInterest, 'YTickLabel', P.sitesOfInterest, 'Box', 'off');
    xlabel(S_fig.hAx, sprintf(xLabel, S_fig.maxAmp));
    ylabel(S_fig.hAx, sprintf(yLabel, S_fig.maxAmp));
    title_(S_fig.hAx, figTitle);
    displayFeature = P.displayFeature;
    S_fig = mergeStructs(S_fig, ...
        makeStruct_(figTitle, primaryCluster, secondaryCluster, sitesOfInterest, xLabel, yLabel, displayFeature));
    S_fig.csHelp = { ...
        '[D]raw polygon', ...
        '[S]plit cluster', ...
        '(shift)+Up/Down: change scale', ...
        '[R]eset scale', ...
        'Zoom: mouse wheel', ...
        'Drag while pressing wheel: pan'};
    set(hFig, 'UserData', S_fig);
end %func
