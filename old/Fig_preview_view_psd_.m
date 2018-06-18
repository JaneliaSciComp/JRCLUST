%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_view_psd_(hFig, vcMode, hMenu)
    % export S_fig.vnThresh and sets vcFile_thresh
    S_fig = get(hFig, 'UserData');
    S_fig.vcPsd_view = vcMode;
    menu_checkbox_(hMenu, vcMode);
    set(hFig, 'UserData', S_fig);
    Fig_preview_plot_([], 1);
end %func
