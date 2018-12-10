function hFigProj = doPlotFigProj(hFigProj, hClust, hCfg)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    iClu1 = S0.iCluCopy;
    iClu2 = S0.iCluPaste;
    update_plot2_proj_(); %erase prev objects

    %---------------
    % Compute
    iSite1 = hClust.viSite_clu(iClu1);
    
    nSitesFigProj = get_set_(hCfg, 'nSitesFigProj', 5); % by request
    nSites = min(nSitesFigProj, size(hCfg.miSites, 1));

    if ~isfield(hCfg, 'viSites_show')
        % P.viSites_show = sort(P.miSites(:, iSite1), 'ascend');

        % center sites around cluster center site
        if nSites < size(hCfg.miSites, 1)
            hCfg.viSites_show = iSite1:iSite1 + nSites - 1;
            if hCfg.viSites_show(end) > max(hCfg.viSite2Chan) % correct for overshooting
                hCfg.viSites_show = hCfg.viSites_show - max(hCfg.viSites_show) + max(hCfg.viSite2Chan);
            end
        else
            hCfg.sitesOfInterest = sort(hCfg.miSites(:, iSite1), 'ascend');
        end
    end

    viSites_show = hCfg.viSites_show;

    cell_plot = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};

    switch lower(hCfg.vcFet_show)
        case {'vpp', 'vmin', 'vmax'}
            vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
            vcYLabel = 'Site # (%0.0f \\muV_{min})';
            
        case {'kilosort', 'pca', 'gpca', 'ppca'}
            S0.pcPair = get_set_(S0, 'pcPair', [1 2]);

            vcXLabel = sprintf('Site # (PC %d)', S0.pcPair(1));
            vcYLabel = sprintf('Site # (PC %d)', S0.pcPair(2));

        otherwise
            vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.vcFet_show, hCfg.vcFet_show, hCfg.vcFet_show);
            vcYLabel = sprintf('Site # (%%0.0f %s)', hCfg.vcFet_show);
    end
    vcTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

    %----------------
    % display
    if isempty(hFigProj.figData)
        hFigProj.figData.maxAmp = hCfg.maxAmp;
        hFigProj.figData.hAx = hFigProj.axes();
        set(hFigProj.figData.hAx, 'Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');
        hFigProj.figData.hPlot0 = line(nan, nan, 'Color', hCfg.mrColor_proj(1,:), 'Parent', hFigProj.figData.hAx);
        hFigProj.figData.hPlot1 = line(nan, nan, 'Color', hCfg.mrColor_proj(2,:), 'Parent', hFigProj.figData.hAx); %place holder
        hFigProj.figData.hPlot2 = line(nan, nan, 'Color', hCfg.mrColor_proj(3,:), 'Parent', hFigProj.figData.hAx); %place holder
        set([hFigProj.figData.hPlot0, hFigProj.figData.hPlot1, hFigProj.figData.hPlot2], cell_plot{:}); %common style
        hFigProj.figData.viSites_show = []; %so that it can update
        hFigProj.figData.vcFet_show = 'vpp';
        % plot boundary
        plotTable_([0, nSites], '-', 'Color', [.5 .5 .5]); %plot in one scoop
        plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5); %plot in one scoop
        mouse_figure(hFigProj);
        set(hFigProj, 'KeyPressFcn', @keyPressFcn_FigProj_);
        hFigProj.figData.hidePlots = mouse_hide_(hFigProj, hFigProj.figData.hPlot0, hFigProj.figData);
        set_fig_(hFigProj, hFigProj.figData);
    end

    % get features for x0,y0,S_plot0 in one go
    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, hCfg.viSites_show);

    if ~isfield(hFigProj.figData, 'viSites_show')
        hFigProj.figData.viSites_show = [];
    end

    if ~equal_vr_(hFigProj.figData.viSites_show, hCfg.viSites_show) || ...
        ~equal_vr_(hFigProj.figData.vcFet_show, hCfg.viSites_show)
        plot_proj_(hFigProj.figData.hPlot0, mrMin0, mrMax0, hCfg, hFigProj.figData.maxAmp);
    end

    plot_proj_(hFigProj.figData.hPlot1, mrMin1, mrMax1, hCfg, hFigProj.figData.maxAmp);
    if ~isempty(iClu2)
        plot_proj_(hFigProj.figData.hPlot2, mrMin2, mrMax2, hCfg, hFigProj.figData.maxAmp);
        vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', iClu1, iClu2, vcTitle);
    else
        update_plot_(hFigProj.figData.hPlot2, nan, nan);
        vcTitle = sprintf('Clu%d (black); %s', iClu1, vcTitle);
    end

    % Annotate axes
    axis_(hFigProj.figData.hAx, [0 nSites 0 nSites]);
    set(hFigProj.figData.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', hCfg.viSites_show, 'YTickLabel', hCfg.viSites_show, 'Box', 'off');
    xlabel(hFigProj.figData.hAx, sprintf(vcXLabel, hFigProj.figData.maxAmp));
    ylabel(hFigProj.figData.hAx, sprintf(vcYLabel, hFigProj.figData.maxAmp));

    % set fig data
    title_(hFigProj.figData.hAx, vcTitle);

    vcFet_show = hCfg.vcFet_show;

    hFigProj.figData = struct_merge_(hFigProj.figData, makeStruct_(vcTitle, iClu1, iClu2, viSites_show, vcXLabel, vcYLabel, vcFet_show));
    hFigProj.figData.csHelp = {'[D]raw polygon', ...
                    '[S]plit cluster', ...
                    '(shift)+Up/Down: change scale', ...
                    '[R]eset scale', ...
                    'Zoom: mouse wheel', ...
                    'Drag while pressing wheel: pan'};
    set(hFigProj, 'UserData', hFigProj.figData);
end %func

%% LOCAL FUNCTIONS
function plot_proj_(hPlot, mrMin, mrMax, P, maxAmp)
    if nargin < 5
        [~, hFigProj.figData] = get_fig_cache_('FigProj');
        maxAmp = hFigProj.figData.maxAmp;
    end

    switch lower(P.vcFet_show)
        case {'vpp', 'vmin', 'vmax'}
            bounds = maxAmp*[0 1];

        otherwise
            % round up to nearest 50 on either side of 0
            bounds = maxAmp*[-1 1];
    end

    [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, bounds, P.maxSite, P);

    % make struct
    maxPair = P.maxSite;
    viSites_show = P.viSites_show;
    S_plot = makeStruct_(mrMax, mrMin, viSites_show, viPlot, tr_dim, maxPair, maxAmp);

    update_plot_(hPlot, vrX, vrY, S_plot);
end %func
