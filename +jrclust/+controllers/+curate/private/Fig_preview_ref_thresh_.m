%--------------------------------------------------------------------------
% 8/17/17 JJJ: created
function Fig_preview_ref_thresh_(hFig)
    S_fig = get(hFig, 'UserData');

    % Ask user and update
    vc_blank_thresh = num2str(S_fig.blankThresh);
    vc_blank_period_ms = num2str(S_fig.blank_period_ms);
    % P = get0_('P');
    % nSites = numel(P.viSite2Chan);
    csAns = inputdlg({'blankThresh (MAD)', 'blank_period_ms (millisecond)'}, ...
    'Common reference threshold', 1, {vc_blank_thresh, vc_blank_period_ms});
    if isempty(csAns), return; end
    % if isnan(site_start) || isnan(site_end), return; end;
    blankThresh = str2num(csAns{1});
    blank_period_ms = str2num(csAns{2});
    if isempty(blankThresh), blankThresh = 0; end
    if isnan(blankThresh) || isnan(blank_period_ms), return; end

    S_fig = set_(S_fig, 'blankThresh', blankThresh, 'blank_period_ms', blank_period_ms);
    set(hFig, 'UserData', S_fig);
    S_fig = doUpdateFigPreview(hFig, S_fig, 1);
end %func
