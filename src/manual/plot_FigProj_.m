%--------------------------------------------------------------------------
function plot_FigProj_(S0)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    S_clu = S0.S_clu;
    P = S0.P;

    [hFig, S_fig] = get_fig_cache_('FigProj');

    iClu1 = S0.iCluCopy;
    iClu2 = S0.iCluPaste;
    update_plot2_proj_(); %erase prev objects

    %---------------
    % Compute
    iSite1 = S_clu.viSite_clu(iClu1);
    
    nSitesFigProj = get_set_(P, 'nSitesFigProj', 5); % by request
    nSites = min(nSitesFigProj, size(P.miSites, 1));

    if ~isfield(P, 'viSites_show')
        % P.viSites_show = sort(P.miSites(:, iSite1), 'ascend');

        % center sites around cluster center site
        if nSites < size(P.miSites, 1)
            P.viSites_show = iSite1:iSite1 + nSites - 1;
            if P.viSites_show(end) > max(P.viSite2Chan) % correct for overshooting
                P.viSites_show = P.viSites_show - max(P.viSites_show) + max(P.viSite2Chan);
            end
        else
            P.sitesOfInterest = sort(P.miSites(:, iSite1), 'ascend');
        end
    end

    viSites_show = P.viSites_show;

    cell_plot = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};

    switch lower(P.vcFet_show)
        case {'vpp', 'vmin', 'vmax'}
            vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
            vcYLabel = 'Site # (%0.0f \\muV_{min})';
            
        case {'kilosort', 'pca', 'gpca', 'ppca'}
            S0.pcPair = get_set_(S0, 'pcPair', [1 2]);

            vcXLabel = sprintf('Site # (PC %d)', S0.pcPair(1));
            vcYLabel = sprintf('Site # (PC %d)', S0.pcPair(2));

        otherwise
            vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.vcFet_show, P.vcFet_show, P.vcFet_show);
            vcYLabel = sprintf('Site # (%%0.0f %s)', P.vcFet_show);
    end
    vcTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

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
        S_fig.viSites_show = []; %so that it can update
        S_fig.vcFet_show = 'vpp';
        % plot boundary
        plotTable_([0, nSites], '-', 'Color', [.5 .5 .5]); %plot in one scoop
        plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5); %plot in one scoop
        mouse_figure(hFig);
        set(hFig, 'KeyPressFcn', @keyPressFcn_FigProj_);
        S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
        set_fig_(hFig, S_fig);
    end

    % get features for x0,y0,S_plot0 in one go
    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, P.viSites_show);

    if ~isfield(S_fig, 'viSites_show')
        S_fig.viSites_show = [];
    end

    if ~equal_vr_(S_fig.viSites_show, P.viSites_show) || ...
        ~equal_vr_(S_fig.vcFet_show, P.viSites_show)
        plot_proj_(S_fig.hPlot0, mrMin0, mrMax0, P, S_fig.maxAmp);
    end

    plot_proj_(S_fig.hPlot1, mrMin1, mrMax1, P, S_fig.maxAmp);
    if ~isempty(iClu2)
        plot_proj_(S_fig.hPlot2, mrMin2, mrMax2, P, S_fig.maxAmp);
        vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', iClu1, iClu2, vcTitle);
    else
        update_plot_(S_fig.hPlot2, nan, nan);
        vcTitle = sprintf('Clu%d (black); %s', iClu1, vcTitle);
    end

    % Annotate axes
    axis_(S_fig.hAx, [0 nSites 0 nSites]);
    set(S_fig.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', P.viSites_show, 'YTickLabel', P.viSites_show, 'Box', 'off');
    xlabel(S_fig.hAx, sprintf(vcXLabel, S_fig.maxAmp));
    ylabel(S_fig.hAx, sprintf(vcYLabel, S_fig.maxAmp));

    % set fig data
    title_(S_fig.hAx, vcTitle);

    vcFet_show = P.vcFet_show;

    S_fig = struct_merge_(S_fig, makeStruct_(vcTitle, iClu1, iClu2, viSites_show, vcXLabel, vcYLabel, vcFet_show));
    S_fig.csHelp = {'[D]raw polygon', ...
                    '[S]plit cluster', ...
                    '(shift)+Up/Down: change scale', ...
                    '[R]eset scale', ...
                    'Zoom: mouse wheel', ...
                    'Drag while pressing wheel: pan'};
    set(hFig, 'UserData', S_fig);
end %func

%% local functions
function plot_proj_(hPlot, mrMin, mrMax, P, maxAmp)
    if nargin < 5
        [~, S_fig] = get_fig_cache_('FigProj');
        maxAmp = S_fig.maxAmp;
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
