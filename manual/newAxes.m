%--------------------------------------------------------------------------
function hAx = newAxes(hFig)
    % clear figure and instantiate a new Axes

    if ischar(hFig)
        hFig = figuresByTag(hFig);
    end

    figure(hFig); % set focus to figure (might be slow)
    clf(hFig);
    hAx = axes(); % create default axes in current figure (hFig)
    hold(hAx, 'on');
end % func
