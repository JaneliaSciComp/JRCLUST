%--------------------------------------------------------------------------
function update_plot2_proj_(newX, newY)
    if nargin == 0
        newX = nan;
        newY = nan;
    end

    [hFig, figData] = getCachedFig('FigProj');

    % erase polygon
    if nargin == 0
        try
            updatePlot(figData.hPlotFG2, newX, newY);
            delete(findobj(get(figData.hAx, 'Child'), 'Type', 'hggroup'));
        catch
            ;
        end
    end
end
