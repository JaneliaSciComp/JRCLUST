%--------------------------------------------------------------------------
function Fig_preview_trange_(hFig, vc_trange, mh)
    % Sets a display time range

    if nargin<1, hFig = []; end
    if nargin<2, vc_trange = 'custom'; end
    if isempty(hFig), hFig = get_fig_cache_('Fig_preview'); end

    if strcmpi(vc_trange, 'custom') % ask user input box
        vcAns = inputdlg('Display time range (s)', 'Time range in seconds', 1, {'.2'});
        if isempty(vcAns), return; end
        trange = str2double(vcAns{1});
    else
        trange = str2double(vc_trange);
    end
    if isnan(trange), return; end
    menu_checkbox_(mh, vc_trange);

    S_fig = get(hFig, 'UserData');
    P = get0_('P');
    S_fig.windowWidth = round(trange * P.sRateHz);
    windowBounds = S_fig.windowBounds(1) + [0, S_fig.windowWidth-1];
    if windowBounds(1)<1
        windowBounds = [1, S_fig.windowWidth];
    elseif windowBounds(2) > S_fig.nSamplesTotal
        windowBounds = [-S_fig.windowWidth+1, 0] + S_fig.nSamplesTotal;
    end
    S_fig.windowBounds = windowBounds;
    set(hFig, 'UserData', S_fig);
    doPlotFigPreview(P);
end %func
