%--------------------------------------------------------------------------
function Fig_preview_site_range_(hFig)
    S_fig = get(hFig, 'UserData');

    % Ask user and update
    vcSiteFrom = num2str(S_fig.siteLim(1));
    vcSiteTo = num2str(S_fig.siteLim(end));
    P = get0_('P');
    nSites = numel(P.viSite2Chan);
    csAns = inputdlg({'Show site from (>=1)', sprintf('Show site to (<=%d)', nSites)}, ...
    'Display site range', 1, {vcSiteFrom, vcSiteTo});
    if isempty(csAns), return; end
    % if isnan(site_start) || isnan(site_end), return; end;
    site_start = max(str2num(csAns{1}), 1);
    site_end = min(str2num(csAns{2}), nSites);
    S_fig = set_(S_fig, 'siteLim', [site_start, site_end]);
    set(S_fig.hAxTraces, 'YLim', S_fig.siteLim + [-1, 1]);
    set(S_fig.hAx_sites, 'YLim', S_fig.siteLim + [-1, 1]);
    set(hFig, 'UserData', S_fig);
end %func
