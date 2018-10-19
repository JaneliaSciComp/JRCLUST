%--------------------------------------------------------------------------
function Fig_preview_spk_thresh_(hFig)
    % update spike threshold qqFactor
    S_fig = get(hFig, 'UserData');
    P = get0_('P');

    % ask user
    vcAns = inputdlg_('Set a spike detection threshold (qqFactor)', 'Spike detection threshold', 1, {num2str(S_fig.qqFactor)});
    if isempty(vcAns), return; end
    qqFactor = str2num(vcAns{1});
    if isnan(qqFactor) || isempty(qqFactor), return; end

    S_fig = Fig_preview_update_(hFig, setfield(S_fig, 'qqFactor', qqFactor), 1);
end %func
