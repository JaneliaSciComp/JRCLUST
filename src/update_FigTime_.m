function update_FigTime_()
    % display features in a new site

    [hFig, S_fig] = get_fig_cache_('FigTime');
    S0 = get(0, 'UserData');
    P = S0.P;
    if ~isVisible_(S_fig.hAx)
        return;
    end

    set0_(P);
    [vrFet0, vrTime0, vcYlabel] = getDispFeaturesSite(S_fig.iSite, [], S0);
    if ~isfield(S_fig, 'fPlot0')
        S_fig.fPlot0 = 1;
    end

    toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);
    update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
    set(S_fig.hPlot1, 'YData', getDispFeaturesSite(S_fig.iSite, S0.iCluCopy, S0));

    if ~isempty(S0.iCluPaste)
        set(S_fig.hPlot2, 'YData', getDispFeaturesSite(S_fig.iSite, S0.iCluPaste, S0));
    else
        hide_plot_(S_fig.hPlot2);
    end

    ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
    imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);
    grid(S_fig.hAx, 'on');
    ylabel(S_fig.hAx, vcYlabel);
end
