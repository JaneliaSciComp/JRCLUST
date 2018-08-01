%--------------------------------------------------------------------------
function traces_test_(P)
    drawnow;
    % csCmd = {'Mouse', 'Menu', 'FigWav', 'FigTime', 'FigWavCor', 'FigProj', 'Exit'};

    % for iCmd = 1:numel(csCmd)
    % vcCmd1 = csCmd{iCmd};
    % fprintf('\tTesting manual-mode %d/%d: %s\n', iCmd, numel(csCmd), vcCmd1);
    hFig = figureByTag('Fig_traces');
    keyPress_fig_(hFig, get_keyPress_('all'));
    try
        close(hFig); %close traces figure. other figures may remain
        close(figureByTag('FigPsd'));
    catch
    end
end %func
