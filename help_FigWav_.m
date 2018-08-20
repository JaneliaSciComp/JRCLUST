%--------------------------------------------------------------------------
function help_FigWav_(hObject, event)
    [~, S_fig] = getCachedFig('FigWav');
    msgbox_(S_fig.csHelp, 1);
end %func
