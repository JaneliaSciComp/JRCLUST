%--------------------------------------------------------------------------
function Fig_preview_trange_(hFig, vc_trange, mh)
    % Sets a display time range

    if nargin<1, hFig = []; end
    if nargin<2, vc_trange = 'custom'; end
    if isempty(hFig), hFig = get_fig_cache_('Fig_preview'); end

    if strcmpi(vc_trange, 'custom') % ask user input box
        vcAns = inputdlg_('Display time range (s)', 'Time range in seconds', 1, {'.2'});
        if isempty(vcAns), return; end
        trange = str2double(vcAns{1});
    else
        trange = str2double(vc_trange);
    end
    if isnan(trange), return; end
    menu_checkbox_(mh, vc_trange);

    S_fig = get(hFig, 'UserData');
    P = get0_('P');
    S_fig.nLoad_bin = round(trange * P.sRateHz);
    nlim_bin = S_fig.nlim_bin(1) + [0, S_fig.nLoad_bin-1];
    if nlim_bin(1)<1
        nlim_bin = [1, S_fig.nLoad_bin];
    elseif nlim_bin(2) > S_fig.nSamples_bin
        nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
    end
    S_fig.nlim_bin = nlim_bin;
    set(hFig, 'UserData', S_fig);
    Fig_preview_plot_(P);
end %func
