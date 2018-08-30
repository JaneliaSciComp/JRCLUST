%--------------------------------------------------------------------------
function S0 = manualMerge(S0)
    if nargin < 1 || isempty(S0)
        S0 = get(0, 'UserData');
    end
    P = S0.P;

    cluster1 = S0.primarySelectedCluster;
    cluster2 = S0.secondarySelectedCluster;

    if isempty(cluster2) || cluster1 == cluster2
        msgbox_('Nothing to merge.', 1);
        return;
    end

    figure_wait_(1);
    drawnow;
    fprintf('%s [W] merging Clu %d and %d\n', datestr(now, 'HH:MM:SS'), cluster1, cluster2);

    % here's where the magic happens
    S_clu = mergePair(S0.S_clu, P, cluster1, cluster2);
    dialogAssert(clusterDataConsistent(S_clu), 'Cluster number is inconsistent after merging');
    S0.S_clu = S_clu;

    set(0, 'UserData', S0);
    plotFigWav(S0); % redraw plot

    S0.primarySelectedCluster = min(cluster1, cluster2);
    S0.secondarySelectedCluster = [];
    set(0, 'UserData', S0);

    updatePlot(S0.hPaste, nan, nan);
    S0 = update_FigCor_(S0);
    S0 = button_CluWav_simulate_(S0.primarySelectedCluster, [], S0);
    S0 = save_log_(sprintf('merge %d %d', S0.primarySelectedCluster, S0.secondarySelectedCluster), S0);
    set(0, 'UserData', S0);

    figure_wait_(0);
    % S_clu = S0.S_clu;
end % function
