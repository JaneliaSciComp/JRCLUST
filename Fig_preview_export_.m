%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function Fig_preview_export_(hFig)
    % Export raw and filtered traces

    S_fig = get(hFig, 'UserData');
    [mnWav_filt, mnWav_raw, S_preview] = deal(S_fig.mnWav_filt, S_fig.mnWav_raw, S_fig.S_preview);
    assignWorkspace_(mnWav_filt, mnWav_raw, S_preview);
    msgbox_('Exported to workspace: Raw traces (mnWav_raw), filtered traces (mnWav_filt), Preview info (S_preview)');
end % function
