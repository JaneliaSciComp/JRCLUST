%--------------------------------------------------------------------------
function updateFigTime()
    % display features in a new site

    [hFig, S_fig] = getCachedFig('FigTime');
    S0 = get(0, 'UserData');
    P = S0.P;
    if ~isVisible_(S_fig.hAx), return ;end
    % P.displayFeature = S_fig.csFet{S_fig.iFet};
    setUserData(P);
    [vrFet0, vrTime0, vcYlabel] = getFet_site_(S_fig.iSite, [], S0);
    if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
    toggleVisible_(S_fig.hPlotBG, S_fig.fPlot0);
    updatePlot(S_fig.hPlotBG, vrTime0, vrFet0);
    set(S_fig.hPlotFG, 'YData', getFet_site_(S_fig.iSite, S0.primarySelectedCluster, S0));
    % set(S_fig.hPlotBG, 'XData', vrTime0, 'YData', vrFet0);
    if ~isempty(S0.secondarySelectedCluster)
        set(S_fig.hPlotFG2, 'YData', getFet_site_(S_fig.iSite, S0.secondarySelectedCluster, S0));
    else
        clearPlots(S_fig.hPlotFG2);
        %     set(S_fig.hPlotFG2, 'XData', nan, 'YData', nan);
    end
    % switch lower(P.displayFeature)
    %     case 'vpp'
    ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
    imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);
    %     otherwise
    %         ylim_(S_fig.hAx, [0, 1] * P.maxAmp);
    %         imrect_set_(S_fig.hRect, [], [0, 1] * P.maxAmp);
    % end
    grid(S_fig.hAx, 'on');
    ylabel(S_fig.hAx, vcYlabel);
end %func
