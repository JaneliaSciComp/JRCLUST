%--------------------------------------------------------------------------
function help_FigWav_(hObject, event)
    [~, S_fig] = get_fig_cache_('FigWav');
    jrclust.utils.qMsgBox(S_fig.helpText, 1);
end %func
