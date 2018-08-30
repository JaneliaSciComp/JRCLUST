%--------------------------------------------------------------------------
function S0 = manualMerge(S0)
    if nargin<1, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;

    if isempty(S0.secondarySelectedCluster)
        msgbox_('Right-click a cluster to merge.', 1); return;
    end
    if S0.primarySelectedCluster == S0.secondarySelectedCluster
        msgbox_('Cannot merge to itself.', 1); return;
    end

    figure_wait_(1); drawnow;
    S0.S_clu = merge_clu_(S0.S_clu, S0.primarySelectedCluster, S0.secondarySelectedCluster, P);
    set(0, 'UserData', S0);
    plotFigWav(S0); %redraw plot
    S0.primarySelectedCluster = min(S0.primarySelectedCluster, S0.secondarySelectedCluster);
    S0.secondarySelectedCluster = [];
    set(0, 'UserData', S0);
    updatePlot(S0.hPaste, nan, nan);
    S0 = update_FigCor_(S0);
    S0 = button_CluWav_simulate_(S0.primarySelectedCluster, [], S0);
    S0 = save_log_(sprintf('merge %d %d', S0.primarySelectedCluster, S0.secondarySelectedCluster), S0);
    set(0, 'UserData', S0);

    % msgbox_close(hMsg);
    figure_wait_(0);
    % S_clu = S0.S_clu;
end % function
