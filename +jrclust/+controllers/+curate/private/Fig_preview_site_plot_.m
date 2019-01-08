%--------------------------------------------------------------------------
function Fig_preview_site_plot_(hFig, vcMode, hMenu)
    S_fig = get(hFig, 'UserData');
    S_fig.siteView = vcMode;
    menu_checkbox_(hMenu, vcMode);
    set(hFig, 'UserData', S_fig);
    doPlotFigPreview([], 1);
end %func
