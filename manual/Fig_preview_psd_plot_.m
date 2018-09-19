%--------------------------------------------------------------------------
% 8/12/17 JJJ: Created
function Fig_preview_psd_plot_(hFig, vcMode, hMenu)
    S_fig = get(hFig, 'UserData');
    switch vcMode
        %     case 'Power'
        %     case 'Detrended'
        case 'Linear', set(S_fig.hAx_psd, 'XScale', 'linear');
        case 'Log', set(S_fig.hAx_psd, 'XScale', 'log');
        otherwise, disperr_(vcMode);
    end %switch
    menu_checkbox_(hMenu, vcMode);
end %func
