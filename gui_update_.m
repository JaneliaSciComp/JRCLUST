%--------------------------------------------------------------------------
function S0 = gui_update_(S0, S_clu)
    if nargin<1, S0 = get(0, 'UserData'); end
    if nargin<2, S_clu = S0.S_clu; end
    plotFigWav(S0); %redraw plot
    S0.primarySelectedCluster = min(S0.primarySelectedCluster, S_clu.nClusters);

    % S0.primarySelectedCluster = 1;
    S0.secondarySelectedCluster = [];
    set(0, 'UserData', S0);
    updatePlot(S0.hPaste, nan, nan); %remove paste cursor
    S0 = update_FigCor_(S0);
    S0 = button_CluWav_simulate_(S0.primarySelectedCluster, [], S0);
    keyPressFcn_cell_(getCachedFig('FigWav'), 'z');
    set(0, 'UserData', S0);
end %func
