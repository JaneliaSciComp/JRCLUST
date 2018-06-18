%--------------------------------------------------------------------------
% deprecated
function plot_SpikePos_(S0, event)


    fPlotX = 0; fPlotAllSites = 0; skip_bg = 4;
    if key_modifier_(event, 'shift'), fPlotX = 1; end
    if key_modifier_(event, 'alt'), fPlotAllSites = 1; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;
    S_clu = S0.S_clu;
    [hFig, S_fig] = get_fig_cache_('FigTime');
    nSites = numel(P.viSite2Chan);
    if nargin<2, fPlotX = 0; end %plot x if set

    if fPlotAllSites
        viSites = 1:nSites;
    else
        iSite_clu = S_clu.viSite_clu(S0.iCluCopy); %only plot spikes from site
        viSites =  P.miSites(:, iSite_clu)';
    end

    % determine spike positinos
    [vrX0, vrY0, vrA0, viClu0, vrT0] = get_spike_clu_(S_clu, viSites); %background
    [vrX2, vrY2, vrA2, viClu2, vrT2] = deal([]);
    if ~fPlotAllSites
        [vrX1, vrY1, vrA1, viClu1, vrT1] = multiindex_(find(viClu0==S0.iCluCopy), ...
        vrX0, vrY0, vrA0, viClu0, vrT0);
        if ~isempty(S0.iCluPaste)
            [vrX2, vrY2, vrA2, viClu2, vrT2] = multiindex_(find(viClu0==S0.iCluPaste), ...
            vrX0, vrY0, vrA0, viClu0, vrT0);
        end
    else
        % display chain
        viClu_Chain(1) = S0.iCluCopy;
        iClu_next = get_next_clu_(S_clu, S0.iCluCopy);
        while ~isempty(iClu_next)
            viClu_Chain(end+1) = iClu_next;
            iClu_next = get_next_clu_(S_clu, iClu_next);
        end
        [vrX1, vrY1, vrA1, viClu1, vrT1] = multiindex_(find(ismember(viClu0, viClu_Chain)), ...
        vrX0, vrY0, vrA0, viClu0, vrT0);
    end

    % if fPlot_ampDist, plot_ampDist_(cmrVpp_site, P); end

    if skip_bg>1 % subsample background
        [vrA0, vrX0, vrY0, vrT0, viClu0] = multiindex_(1:skip_bg:numel(vrA0), ...
        vrA0, vrX0, vrY0, vrT0, viClu0);
    end

    S_fig.xylim = [get(S_fig.hAx, 'XLim'), get(S_fig.hAx, 'YLim')]; %store limits
    S_fig.xylim_track = S_fig.xylim;

    %------------
    % Draw
    if ~isfield(S_fig, 'vhAx_track')
        S_fig.vhAx_track = axes('Parent', hFig, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');
        hold(S_fig.vhAx_track, 'on');
        xlabel(S_fig.vhAx_track, 'Time (s)');
        S_fig.hPlot0_track = scatter(nan, nan, 5, nan, 'filled'); %place holder
        S_fig.hPlot1_track = line(nan, nan, 'Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 5, 'LineStyle', 'none');
        S_fig.hPlot2_track = line(nan, nan, 'Color', [1 0 0], 'Marker', 'o', 'MarkerSize', 5, 'LineStyle', 'none');
    else
        toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 1);
        set(S_fig.vhAx_track, 'Visible', 'on');
    end
    % axes(S_fig.vhAx_track);
    mouse_figure(hFig, S_fig.vhAx_track);
    toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2}, 0);
    if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
    toggleVisible_(S_fig.hPlot0_track, S_fig.fPlot0);
    if fPlotX
        set(S_fig.hPlot0_track, 'XData', vrT0, 'YData', vrX0, 'CData', log10(vrA0));
        update_plot_(S_fig.hPlot1_track, vrT1, vrX1);
        update_plot_(S_fig.hPlot2_track, vrT2, vrX2);
        ylabel(S_fig.vhAx_track, 'X Pos [pix]');
        S_fig.xylim_track(3:4) = [0 1.5];
    else
        set(S_fig.hPlot0_track, 'XData', vrT0, 'YData', vrY0, 'CData', log10(vrA0));
        update_plot_(S_fig.hPlot1_track, vrT1, vrY1);
        update_plot_(S_fig.hPlot2_track, vrT2, vrY2);
        ylabel(S_fig.vhAx_track, 'Y Pos [pix]');
        S_fig.xylim_track(3:4) = round(median(vrY1)) + [-1,1] * floor(P.maxSite);
    end
    axis_(S_fig.vhAx_track, S_fig.xylim_track);
    colormap(S_fig.vhAx_track, flipud(colormap('gray')));
    set(S_fig.vhAx_track, 'CLim', [.5 3.5]);
    grid(S_fig.vhAx_track, 'on');

    % Set title
    if isempty(S0.iCluPaste)
        vcTitle = sprintf('Color: log10 Vpp [uV]; Clu%d(black); Press [T] to return; [B]ackground', S0.iCluCopy);
    else
        vcTitle = sprintf('Color: log10 Vpp [uV]; Clu%d(black); Clu%d(red); Press [T] to return; [B]ackground', S0.iCluCopy, S0.iCluPaste);
    end
    title(S_fig.vhAx_track, vcTitle);
    set(hFig, 'UserData', S_fig);
end %func
