%--------------------------------------------------------------------------
function Fig_preview_site_thresh_(hFig)
    % Site preview threshold
    S_fig = get(hFig, 'UserData');
    P = get0_('P');

    % ask user
    if isempty(S_fig.thresh_corr_bad_site)
        vc_thresh_site = '0';
    else
        vc_thresh_site = num2str(S_fig.thresh_corr_bad_site);
    end
    csAns = inputdlg_('Set a correlation threshold [0-1) for detecting bad sites (0 to disable)', 'thresh_corr_bad_site (set to 0 to disable)', 1, {vc_thresh_site});
    if isempty(csAns), return; end
    thresh_corr_bad_site = str2num(csAns{1});
    if thresh_corr_bad_site>=1 || thresh_corr_bad_site<0 || isnan(thresh_corr_bad_site)
        msgbox_('Invalid range');
        return;
    end

    S_fig.thresh_corr_bad_site = thresh_corr_bad_site;
    set(hFig, 'UserData', S_fig);
    Fig_preview_update_(hFig, S_fig, 1);

end % function
