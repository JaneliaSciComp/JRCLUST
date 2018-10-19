%--------------------------------------------------------------------------
% 9/26/17 JJJ: Bugfix: returning S0 so that the cluster can be updated
function S0 = ui_delete_(S0)
    if nargin<1, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;
    if ~isempty(S0.iCluPaste)
        msgbox_('Must select one cluster', 1); return;
    end
    figure_wait_(1);

    iClu_del = S0.iCluCopy;
    % hMsg = msgbox_open_('Deleting...');
    S0.S_clu = delete_clu_(S0.S_clu, S0.iCluCopy);
    set(0, 'UserData', S0);
    plot_FigWav_(S0); %redraw plot
    % S0.S_clu.mrWavCor = wavCor_delete_(S0.iCluCopy);
    FigWavCor_update_(S0);
    S0.iCluCopy = min(S0.iCluCopy, S0.S_clu.nClu);
    % set(0, 'UserData', S0);
    button_CluWav_simulate_(S0.iCluCopy);

    % close_(hMsg);
    figure_wait_(0);
    fprintf('%s [W] deleted Clu %d\n', datestr(now, 'HH:MM:SS'), iClu_del);
    S0 = save_log_(sprintf('delete %d', iClu_del), S0);
    set(0, 'UserData', S0);
end
