%--------------------------------------------------------------------------
function traces_test_(P)
    drawnow;
    % csCmd = {'Mouse', 'Menu', 'FigWav', 'FigTime', 'FigWavCor', 'FigProj', 'Exit'};

    % for iCmd = 1:numel(csCmd)
    % vcCmd1 = csCmd{iCmd};
    % fprintf('\tTesting manual-mode %d/%d: %s\n', iCmd, numel(csCmd), vcCmd1);
    hFig = get_fig_('Fig_traces');
    keyPress_fig_(hFig, get_keyPress_('all'));
    try
        jrclust.utils.tryClose(hFig); %close traces figure. other figures may remain
        jrclust.utils.tryClose(get_fig_('FigPsd'));
    catch
    end
end %func
