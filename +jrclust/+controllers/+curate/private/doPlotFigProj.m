function hFigProj = doPlotFigProj(hFigProj, hClust, hCfg, selected)
    %DOPLOTFIGPROJ Plot feature projection figure
    update_plot2_proj_(); %erase prev objects

    iSite = hClust.clusterSites(selected(1));
    
    % limit the number of sites to display in the feature projection view
    nSites = min(hCfg.nSitesFigProj, size(hCfg.siteNeighbors, 1)); % by request
    % center sites around cluster center site
    if nSites < size(hCfg.siteNeighbors, 1)
        sitesToShow = iSite:iSite + nSites - 1;
        if sitesToShow(end) > max(hCfg.siteMap) % correct for overshooting
            sitesToShow = sitesToShow - max(sitesToShow) + max(hCfg.siteMap);
        end
    else
        sitesToShow = sort(hCfg.siteNeighbors(:, iSite), 'ascend');
    end

    if strcmp(hCfg.dispFeature, 'vpp')
        xLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        yLabel = 'Site # (%0.0f \\muV_{min})';    
    elseif ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        xLabel = sprintf('Site # (PC %d)', 1); % fix this in post
        yLabel = sprintf('Site # (PC %d)', 2);
    else
        xLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        yLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end
    figTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

    if isempty(hFigProj.figData)
        hFigProj.figData.maxAmp = hCfg.maxAmp;
        hFigProj.axes();
        hFigProj.axSet('Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');

        hFigProj.addLine('background',  nan, nan, 'Color', hCfg.mrColor_proj(1, :));
        hFigProj.addLine('foreground',  nan, nan, 'Color', hCfg.mrColor_proj(2, :)); % placeholder
        hFigProj.addLine('foreground2', nan, nan, 'Color', hCfg.mrColor_proj(3, :)); % placeholder
        
        plotStyle = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};
        hFigProj.plotSet('background', plotStyle);
        hFigProj.plotSet('foreground', plotStyle);
        hFigProj.plotSet('foreground2', plotStyle);

%         hFigProj.figData.sitesToShow = []; %so that it can update
%         hFigProj.figData.dispFeature = 'vpp';

        %%%% PICK UP HERE
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
        ~equal_vr_(hFigProj.figData.dispFeature, hCfg.viSites_show)
        plot_proj_(hFigProj.figData.hPlot0, mrMin0, mrMax0, hCfg, hFigProj.figData.maxAmp);
    end

    plot_proj_(hFigProj.figData.hPlot1, mrMin1, mrMax1, hCfg, hFigProj.figData.maxAmp);
    if numel(selected) == 2
        plot_proj_(hFigProj.figData.hPlot2, mrMin2, mrMax2, hCfg, hFigProj.figData.maxAmp);
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', selected(1), selected(2), figTitle);
    else
        update_plot_(hFigProj.figData.hPlot2, nan, nan);
        figTitle = sprintf('Clu%d (black); %s', selected(1), figTitle);
    end

    % Annotate axes
    axis_(hFigProj.figData.hAx, [0 nSites 0 nSites]);
    set(hFigProj.figData.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', hCfg.viSites_show, 'YTickLabel', hCfg.viSites_show, 'Box', 'off');
    xlabel(hFigProj.figData.hAx, sprintf(xLabel, hFigProj.figData.maxAmp));
    ylabel(hFigProj.figData.hAx, sprintf(yLabel, hFigProj.figData.maxAmp));

    % set fig data
    title_(hFigProj.figData.hAx, figTitle);

    dispFeature = hCfg.dispFeature;

    hFigProj.figData = struct_merge_(hFigProj.figData, makeStruct_(figTitle, selected(1), selected(2), sitesToShow, xLabel, yLabel, dispFeature));
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

    switch lower(P.dispFeature)
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
