%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function Fig_preview_export_(hFig)
    % Export raw and filtered traces

    S_fig = get(hFig, 'UserData');
    [mnWav_filt, tracesRaw, S_preview] = deal(S_fig.mnWav_filt, S_fig.tracesRaw, S_fig.S_preview);
    assignWorkspace_(mnWav_filt, tracesRaw, S_preview);
    jrclust.utils.qMsgBox('Exported to workspace: Raw traces (tracesRaw), filtered traces (mnWav_filt), Preview info (S_preview)');
end %func
