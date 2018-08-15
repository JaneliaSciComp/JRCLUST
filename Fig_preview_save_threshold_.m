%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_save_threshold_(hFig)
    % export S_fig.vnThresh and sets thresholdFile
    P = get0_('P');
    S_fig = get(hFig, 'UserData');
    siteThresholds = S_fig.siteThresholds;
    thresholdFile = strrep(P.paramFile, '.prm', '_thresh.mat');
    save(thresholdFile, 'siteThresholds'); % also need to store filter values?
    P.thresholdFile = thresholdFile;
    setUserData(P);
    updateParamFile(P, P.paramFile);

    msgbox_(sprintf('Saved to %s and updated %s (thresholdFile)', ...
    thresholdFile, P.paramFile));
end %func
