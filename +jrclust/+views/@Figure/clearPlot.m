function clearPlot(obj, plotKey)
    %CLEARPLOT Set XData and YData of a plot to nan
    obj.updatePlot(plotKey, nan, nan); % updatePlot checks for existence of plotKey
end