%--------------------------------------------------------------------------
function help_FigWav_(hObject, event)
    [~, S_fig] = get_fig_cache_('FigWav');
    msgbox_(S_fig.csHelp, 1);
end %func
