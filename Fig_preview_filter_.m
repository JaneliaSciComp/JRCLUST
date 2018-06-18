%--------------------------------------------------------------------------
function Fig_preview_filter_(hFig, vcMode, hMenu)
    % Set reference types
    S_fig = get(hFig, 'UserData');
    menu_checkbox_(hMenu, vcMode);
    S_fig = Fig_preview_update_(hFig, set_(S_fig, 'vcFilter', vcMode, 'fFilter', 1), 1);
end %func
