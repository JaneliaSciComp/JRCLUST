%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_view_ref_(hFig, vcMode, hMenu)
    S_fig = get(hFig, 'UserData');
    S_fig.refView = vcMode;
    menu_checkbox_(hMenu, vcMode);
    set(hFig, 'UserData', S_fig);
    doPlotFigPreview([], 1);
end %func
