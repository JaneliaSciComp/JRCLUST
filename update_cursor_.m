%--------------------------------------------------------------------------
function  S0 = update_cursor_(S0, iClu, fPaste)
    if isempty(iClu), return; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    % [hFig, S_fig] = getCachedFig('FigWav');

    if ~isfield(S0, 'hCopy'), S0.hCopy = []; end
    if ~isfield(S0, 'hPaste'), S0.hPaste = []; end

    if ~fPaste
        primarySelectedCluster = iClu;
        if primarySelectedCluster <1 || primarySelectedCluster > S_clu.nClusters, return; end
        updatePlot(S0.hPaste, nan, nan); %hide paste
        S0.secondarySelectedCluster = [];
        [S0.primarySelectedCluster, S0.hCopy] = plot_tmrWav_clu_(S0, primarySelectedCluster, S0.hCopy, [0 0 0]);
    else
        secondarySelectedCluster = iClu;
        if secondarySelectedCluster < 1 || secondarySelectedCluster > S_clu.nClusters || S0.primarySelectedCluster == secondarySelectedCluster, return; end
        [S0.secondarySelectedCluster, S0.hPaste] = plot_tmrWav_clu_(S0, secondarySelectedCluster, S0.hPaste, [1 0 0]);
    end
    % set(hFig, 'UserData', S_fig);
    cursorFigClusterCor(S0);
    if nargout==0, set(0, 'UserData', S0); end
end %func
