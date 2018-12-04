%--------------------------------------------------------------------------
function plot_FigIsi_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    [hFig, S_fig] = get_fig_cache_('FigIsi');

    [vrX1, vrY1] = get_returnMap_(S0.iCluCopy, P);
    if isempty(S_fig)
        S_fig.hAx = hFig.axes();
        S_fig.hPlot1 = plot(S_fig.hAx, nan, nan, 'ko');
        S_fig.hPlot2 = plot(S_fig.hAx, nan, nan, 'ro');
        set(S_fig.hAx, 'XScale','log', 'YScale','log');
        xlabel('ISI_{k} (ms)'); ylabel('ISI_{k+1} (ms)');
        axis_(S_fig.hAx, [1 10000 1 10000]);
        grid(S_fig.hAx, 'on');
        % show refractory line
        line(get(S_fig.hAx,'XLim'), P.spkRefrac_ms*[1 1], 'Color', [1 0 0]);
        line(P.spkRefrac_ms*[1 1], get(S_fig.hAx,'YLim'), 'Color', [1 0 0]);
    end
    update_plot_(S_fig.hPlot1, vrX1, vrY1);
    if ~isempty(S0.iCluPaste)
        [vrX2, vrY2] = get_returnMap_(S0.iCluPaste, P);
        update_plot_(S_fig.hPlot2, vrX2, vrY2);
    else
        update_plot_(S_fig.hPlot2, nan, nan);
    end

    set(hFig, 'UserData', S_fig);
end %func
