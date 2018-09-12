%--------------------------------------------------------------------------
function ui_update_log_(cS_log, S0)
    % the last one is selected
    % persistent mh_history

    % List recent activities
    % if nargin<2, S0 = get(0, 'UserData'); end
    % set(hMenu_history, 'Label', sprintf('Undo %s', cS_log.csCmd{end}), 'Enable', 'on');
    if nargin<2, S0=[]; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;
    mh_history = get_tag_('mh_history', 'uimenu');

    % Delete children and update
    delete_(mh_history.Children); %kill all children
    for iMenu = 1:numel(cS_log) % reverse order
        iLog = numel(cS_log) - iMenu + 1;
        S_log1 = cS_log{iLog};
        vcLabel1 = sprintf('%s: %s', datestr(S_log1.datenum), S_log1.vcCmd);
        fEnable = (iMenu <= P.MAX_LOG) && iMenu~=1;
        uimenu(mh_history, 'Label', vcLabel1, 'Callback', @(h,e)restore_log_(iMenu), ...
        'Checked', ifeq_(iMenu==1, 'on', 'off'), ...
        'Enable', ifeq_(fEnable, 'on', 'off'));
    end
    % update undo/redo menu
end %func
