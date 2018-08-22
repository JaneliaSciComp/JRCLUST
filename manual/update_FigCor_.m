%--------------------------------------------------------------------------
function S0 = update_FigCor_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    [hFig, S_fig] = getCachedFig('FigClusterCor');

    % figure(hFig);
    % hMsg = msgbox_open('Computing Correlation');

    xylim = get(gca, {'XLim', 'YLim'});
    S0.S_clu = S_clu;
    plotFigClusterCor(S0);
    set(gca, {'XLim', 'YLim'}, xylim);

    set(S_fig.hCursorV, 'XData', S0.primarySelectedCluster*[1 1]);
    set(S_fig.hCursorH, 'YData', S0.primarySelectedCluster*[1 1]);

    if nargout==0
        set(0, 'UserData', S0); %update field
    end
end %func
