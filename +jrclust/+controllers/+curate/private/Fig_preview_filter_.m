%--------------------------------------------------------------------------
function Fig_preview_filter_(hFig, vcMode, hMenu)
    % Set reference types
    S_fig = get(hFig, 'UserData');
    menu_checkbox_(hMenu, vcMode);
    S_fig = doUpdateFigPreview(hFig, set_(S_fig, 'filterType', vcMode, 'fFilter', 1), 1);
end %func
