%--------------------------------------------------------------------------
function S0 = ui_merge_(S0)
    if nargin<1, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;

    if isempty(S0.iCluPaste)
        msgbox_('Right-click a cluster to merge.', 1); return;
    end
    if S0.iCluCopy == S0.iCluPaste
        msgbox_('Cannot merge to itself.', 1); return;
    end

    figure_wait_(1); drawnow;
    S0.S_clu = merge_clu_(S0.S_clu, S0.iCluCopy, S0.iCluPaste, P);
    set(0, 'UserData', S0);
    plot_FigWav_(S0); %redraw plot
    S0.iCluCopy = min(S0.iCluCopy, S0.iCluPaste);
    S0.iCluPaste = [];
    set(0, 'UserData', S0);
    update_plot_(S0.hPaste, nan, nan);
    S0 = update_FigCor_(S0);
    S0 = button_CluWav_simulate_(S0.iCluCopy, [], S0);
    S0 = save_log_(sprintf('merge %d %d', S0.iCluCopy, S0.iCluPaste), S0);
    set(0, 'UserData', S0);

    % msgbox_close(hMsg);
    figure_wait_(0);
    % S_clu = S0.S_clu;
end %func
