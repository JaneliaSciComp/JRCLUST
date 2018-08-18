%--------------------------------------------------------------------------
function  S0 = update_cursor_(S0, iClu, fPaste)
    if isempty(iClu), return; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    % [hFig, S_fig] = getCachedFig('FigWav');

    if ~isfield(S0, 'hCopy'), S0.hCopy = []; end
    if ~isfield(S0, 'hPaste'), S0.hPaste = []; end

    if ~fPaste
        iCluCopy = iClu;
        if iCluCopy <1 || iCluCopy > S_clu.nClusters, return; end
        update_plot_(S0.hPaste, nan, nan); %hide paste
        S0.iCluPaste = [];
        [S0.iCluCopy, S0.hCopy] = plot_tmrWav_clu_(S0, iCluCopy, S0.hCopy, [0 0 0]);
    else
        iCluPaste = iClu;
        if iCluPaste < 1 || iCluPaste > S_clu.nClusters || S0.iCluCopy == iCluPaste, return; end
        [S0.iCluPaste, S0.hPaste] = plot_tmrWav_clu_(S0, iCluPaste, S0.hPaste, [1 0 0]);
    end
    % set(hFig, 'UserData', S_fig);
    cursor_FigClusterCor_(S0);
    if nargout==0, set(0, 'UserData', S0); end
end %func
