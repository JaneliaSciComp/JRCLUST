%--------------------------------------------------------------------------
function restore_log_(iMenu1)
    % persistent mh_history
    figure_wait_(1);
    [cS_log, miClu_log, P] = get0_('cS_log', 'miClu_log', 'P');
    S_clu1 = cS_log{end - iMenu1 + 1}; % last ones shown first
    S_clu1.viClu = int32(miClu_log(:,iMenu1));

    hMsg = msgbox_(sprintf('Restoring to %s (%s)', S_clu1.vcCmd, datestr(S_clu1.datenum)), 0);
    [S_clu, S0] = S_clu_new_(S_clu1);

    % Update checkbox
    mh_history = get_tag_('mh_history', 'uimenu');
    vhMenu = mh_history.Children;
    vhMenu = vhMenu(end:-1:1); %reverse order
    for iMenu = 1:numel(vhMenu)
        fTargetItem = iMenu==iMenu1;
        fEnable = ~fTargetItem && iMenu <= P.MAX_LOG;
        set(vhMenu(iMenu), ...
        'Checked', ifeq_(fTargetItem, 'on', 'off'), ...
        'Enable', ifeq_(fEnable, 'on', 'off'));
    end %for

    % update GUI
    S0 = gui_update_(S0, S_clu);
    % plotFigWav(S0); %redraw plot
    % S0.primarySelectedCluster = min(S0.primarySelectedCluster, S_clu.nClusters);
    % S0.secondarySelectedCluster = [];
    % set(0, 'UserData', S0);
    % updatePlot(S0.hPaste, nan, nan); %remove paste cursor
    % S0 = update_FigCor_(S0);
    % S0 = button_CluWav_simulate_(S0.primarySelectedCluster, [], S0);
    % set(0, 'UserData', S0);
    tryClose(hMsg);
    figure_wait_(0);
end % function
