%--------------------------------------------------------------------------
function Fig_preview_ref_(hFig, vcMode, hMenu)
    % Set reference types
    S_fig = get(hFig, 'UserData');
    S_fig.vcCommonRef = vcMode;
    S_fig = Fig_preview_update_(hFig, S_fig, 1);
    menu_checkbox_(hMenu, vcMode);
end %func
