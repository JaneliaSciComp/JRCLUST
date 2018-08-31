%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_view_ref_(hFig, vcMode, hMenu)
    S_fig = get(hFig, 'UserData');
    S_fig.vcRef_view = vcMode;
    menu_checkbox_(hMenu, vcMode);
    set(hFig, 'UserData', S_fig);
    Fig_preview_plot_([], 1);
end % function
