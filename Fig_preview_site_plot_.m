%--------------------------------------------------------------------------
function Fig_preview_site_plot_(hFig, vcMode, hMenu)
    S_fig = get(hFig, 'UserData');
    S_fig.vcSite_view = vcMode;
    menu_checkbox_(hMenu, vcMode);
    set(hFig, 'UserData', S_fig);
    Fig_preview_plot_([], 1);
end %func
