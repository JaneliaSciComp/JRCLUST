%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function Fig_preview_save_prm_(hFig)
    % Update to a parameter file, show a preview window
    % List of parameters to update. Set threshold: save to a known location and load
    P = get0_('P');
    S_fig = get(hFig, 'UserData');

    % Select variables to export
    [fft_thresh, qqFactor, vlSite_bad, vcFilter, vcCommonRef, blank_thresh, blank_period_ms] = ...
    get_(S_fig, 'fft_thresh', 'qqFactor', 'vlSite_bad', 'vcFilter', 'vcCommonRef', 'blank_thresh', 'blank_period_ms');
    viSiteZero = find(vlSite_bad);
    P_update = makeStruct_(fft_thresh, viSiteZero, qqFactor, vcFilter, vcCommonRef, blank_thresh, blank_period_ms);

    % Preview variables in the edit box
    vcUpdate = struct2str_(P_update);
    csAns = inputdlg_(P.paramFile, 'Update confirmation', 16, {vcUpdate}, struct('Resize', 'on'));
    if isempty(csAns), return; end
    [P_update, vcErr] = str2struct_(csAns{1});
    if isempty(P_update)
        msgbox_(vcErr);
        return;
    end
    P = get0_('P');
    updateParamFile(P_update, P.paramFile);
    % P = mergeStructs(P, P_update);
    % setUserData(P);
    edit(P.paramFile);
end %func
